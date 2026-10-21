"""Find homologous sequences of the target protein with remote search services."""

# The homologs are returned as rows aligned to the target sequence: one column
# per target residue, "-" where a homolog has no aligned residue. The ColabFold
# MMseqs2 server is tried first, the EBI BLAST service (UniProtKB) second.
# Near-identical hits carry no conservation signal (and make WHISCY fail when
# all of them are identical), so both are filtered for redundancy.

import io
import json
import logging
import os
import re
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request

logger = logging.getLogger("whiscy_setup")

COLABFOLD_URL = "https://api.colabfold.com"
EBI_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
USER_AGENT = "whiscy (https://github.com/haddocking/whiscy)"

# Redundancy filter, validated in haddocking/whiscy#39
MAX_IDENTITY = 0.9
MIN_COVERAGE = 0.5
MAX_HITS = 250
# Fewer hits than this after filtering and the next service is tried
MIN_HITS = 10
# Fewer hits than this and WHISCY cannot score the alignment (1 and 2 hits
# failed in the validation)
MIN_USABLE_HITS = 3
# Hits read from each service before filtering
COLABFOLD_CANDIDATES = 5000
EBI_CANDIDATES = 1000

# Time limit for one service, in seconds
SEARCH_TIMEOUT = int(os.environ.get("WHISCY_SEARCH_TIMEOUT", "3600"))
POLL_INTERVAL = 10
HTTP_TIMEOUT = 120
HTTP_TRIES = 5
HTTP_RETRY_WAIT = 30


class HomologSearchError(Exception):
    """No usable homologs were found."""


def _http(url, data=None):
    """GET (or POST data) an https URL, retrying on network errors."""
    if not url.startswith("https://"):
        raise HomologSearchError("Not an https URL: {0}".format(url))
    if data is not None:
        data = urllib.parse.urlencode(data).encode()
    request = urllib.request.Request(
        url, data=data, headers={"User-Agent": USER_AGENT}
    )
    for attempt in range(1, HTTP_TRIES + 1):
        try:
            # Only the fixed https service URLs above are requested
            with urllib.request.urlopen(  # nosec B310 # nosemgrep
                request, timeout=HTTP_TIMEOUT
            ) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if attempt == HTTP_TRIES:
                raise HomologSearchError("{0}: {1}".format(url, err)) from err
            logger.warning("Request failed ({0}), retrying".format(err))
            time.sleep(HTTP_RETRY_WAIT)


def _wait(check, service):
    """Poll check() until it returns a final status or the time limit passes."""
    start = time.time()
    while True:
        status = check()
        if status is not None:
            return status
        if time.time() - start > SEARCH_TIMEOUT:
            raise HomologSearchError(
                "{0} search did not finish within {1} s".format(
                    service, SEARCH_TIMEOUT
                )
            )
        time.sleep(POLL_INTERVAL)


def colabfold_search(sequence):
    """Search UniRef with the ColabFold MMseqs2 server and return the A3M text."""
    query = {"q": ">101\n{0}".format(sequence), "mode": "all"}
    ticket = json.loads(_http("{0}/ticket/msa".format(COLABFOLD_URL), query))

    def check():
        nonlocal ticket
        if ticket["status"] in ("UNKNOWN", "PENDING", "RUNNING"):
            ticket = json.loads(
                _http("{0}/ticket/{1}".format(COLABFOLD_URL, ticket["id"]))
            )
            return None
        return ticket["status"]

    status = _wait(check, "ColabFold")
    if status != "COMPLETE":
        raise HomologSearchError(
            "ColabFold search ended with status {0}".format(status)
        )
    archive = _http("{0}/result/download/{1}".format(COLABFOLD_URL, ticket["id"]))
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        return tar.extractfile("uniref.a3m").read().decode()


def ebi_search(sequence, email):
    """Run blastp against UniProtKB at the EBI and return the result JSON."""
    params = {
        "email": email,
        "program": "blastp",
        "stype": "protein",
        "database": "uniprotkb",
        "sequence": sequence,
        "alignments": EBI_CANDIDATES,
        "scores": EBI_CANDIDATES,
        "exp": "10",
    }
    job = _http("{0}/run".format(EBI_URL), params).decode()

    def check():
        status = _http("{0}/status/{1}".format(EBI_URL, job)).decode()
        return None if status in ("QUEUED", "RUNNING") else status

    status = _wait(check, "EBI BLAST")
    if status != "FINISHED":
        raise HomologSearchError(
            "EBI BLAST search ended with status {0}".format(status)
        )
    return _http("{0}/result/{1}/json".format(EBI_URL, job)).decode()


def parse_a3m(text, sequence):
    """Rows of an A3M alignment (query first) aligned to the target sequence."""
    rows = []
    name = None
    chunks = []
    # The server separates the alignments of several queries with null bytes
    for line in text.replace("\x00", "").splitlines() + [">"]:
        if line.startswith(">"):
            if name is not None:
                # Lower case letters are insertions relative to the query
                rows.append((name, re.sub("[a-z.]", "", "".join(chunks))))
            name = line[1:].split()[0] if len(line) > 1 else None
            chunks = []
        elif line.strip() and not line.startswith("#"):
            chunks.append(line.strip())
    hits = rows[1:]
    return [(name, row) for name, row in hits if len(row) == len(sequence)]


def parse_ebi_json(text, sequence):
    """Hits of an EBI BLAST JSON result aligned to the target sequence."""
    rows = []
    for hit in json.loads(text).get("hits", []):
        hsps = hit.get("hit_hsps") or []
        if not hsps:
            continue
        hsp = hsps[0]
        row = anchor_on_query(
            hsp["hsp_qseq"], hsp["hsp_hseq"], int(hsp["hsp_query_from"]), sequence
        )
        rows.append((hit.get("hit_acc") or hit.get("hit_id"), row))
    return rows


def anchor_on_query(query_aligned, hit_aligned, query_start, sequence):
    """Place a pairwise alignment on the target's columns (query gaps dropped)."""
    row = ["-"] * len(sequence)
    position = query_start - 1
    for q, h in zip(query_aligned, hit_aligned):
        if q == "-":
            continue
        if 0 <= position < len(sequence):
            row[position] = h
        position += 1
    return "".join(row)


def _identity(a, b):
    """Fraction of identical residues over the columns where both have one."""
    pairs = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    if not pairs:
        return 1.0
    return sum(x == y for x, y in pairs) / len(pairs)


def filter_redundant(
    sequence,
    rows,
    max_identity=MAX_IDENTITY,
    min_coverage=MIN_COVERAGE,
    max_hits=MAX_HITS,
):
    """Keep, in order, the hits below max_identity to the target and to each other."""
    kept = []
    for name, row in rows:
        if sum(c != "-" for c in row) < min_coverage * len(sequence):
            continue
        if _identity(sequence, row) > max_identity:
            continue
        if any(_identity(other, row) > max_identity for _, other in kept):
            continue
        kept.append((name, row))
        if len(kept) == max_hits:
            break
    return kept


def find_homologs(sequence, email=None):
    """Diverse homologs of the sequence: ColabFold first, then the EBI."""
    services = [
        (
            "ColabFold MMseqs2",
            lambda: parse_a3m(colabfold_search(sequence), sequence),
        )
    ]
    if email:
        services.append(
            (
                "EBI BLAST",
                lambda: parse_ebi_json(ebi_search(sequence, email), sequence),
            )
        )
    else:
        logger.warning("No contact email set, the EBI BLAST fallback is disabled")

    best = []
    for service, search in services:
        logger.info("Searching homologs with {0}".format(service))
        try:
            hits = search()[:COLABFOLD_CANDIDATES]
        except (
            HomologSearchError,
            AttributeError,
            KeyError,
            TypeError,
            ValueError,
            tarfile.TarError,
        ) as err:
            logger.warning("{0} failed: {1}".format(service, err))
            continue
        kept = filter_redundant(sequence, hits)
        logger.info(
            "{0}: {1} hits, {2} kept (redundant ones removed)".format(
                service, len(hits), len(kept)
            )
        )
        if len(kept) >= MIN_HITS:
            return kept
        if len(kept) > len(best):
            best = kept
    if len(best) < MIN_USABLE_HITS:
        raise HomologSearchError(
            "Not enough homologs found for this sequence ({0}). "
            "Please upload a multiple sequence alignment.".format(len(best))
        )
    logger.warning(
        "Only {0} homologs found, predictions will be less reliable".format(len(best))
    )
    return best
