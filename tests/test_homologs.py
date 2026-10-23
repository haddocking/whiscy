"""Tests for the remote homolog search."""

import io
import json
import tarfile
import urllib.error

import pytest

from whiscy.modules import homologs

QUERY = "MKTAYIAKQR"

A3M = """#10\t1
>101
MKTAYIAKQR
>UniRef100_A
MKTAYIAKQK
>UniRef100_B
MKaTAYIGKQR
>UniRef100_C
--TAYLAKQR
\x00"""

EBI_JSON = json.dumps(
    {
        "hits": [
            {
                "hit_acc": "P1",
                "hit_id": "P1_HUMAN",
                "hit_hsps": [
                    {
                        "hsp_query_from": 3,
                        "hsp_qseq": "TAY-IAKQR",
                        "hsp_hseq": "TAYWIAKQK",
                    }
                ],
            },
            {
                "hit_acc": "Q2",
                "hit_hsps": [
                    {
                        "hsp_query_from": 1,
                        "hsp_qseq": "MKTAYIAKQR",
                        "hsp_hseq": "MKT-YIGKQR",
                    }
                ],
            },
            {"hit_acc": "R3", "hit_hsps": []},
        ]
    }
)


def test_parse_a3m():
    """Insertions are removed and the query row is skipped."""
    rows = homologs.parse_a3m(A3M, QUERY)
    assert rows == [
        ("UniRef100_A", "MKTAYIAKQK"),
        ("UniRef100_B", "MKTAYIGKQR"),
        ("UniRef100_C", "--TAYLAKQR"),
    ]


def test_parse_ebi_json():
    """Hits are placed on the query columns; query gaps are dropped."""
    rows = homologs.parse_ebi_json(EBI_JSON, QUERY)
    assert rows == [("P1", "--TAYIAKQK"), ("Q2", "MKT-YIGKQR")]


def test_anchor_on_query_clips_to_the_sequence():
    """Columns beyond the target sequence are ignored."""
    assert homologs.anchor_on_query("KQRXX", "KQKAA", 8, QUERY) == "-------KQK"


def test_filter_redundant():
    """Near-identical, redundant and short hits are removed, order is kept."""
    rows = [
        ("same", QUERY),  # 100% identical to the query
        ("a", "MKTSYLAKER"),  # 70%
        ("a_copy", "MKTSYLAKER"),  # identical to a
        ("short", "MKTS------"),  # covers 40% of the query
        ("b", "LRTAWIGKQE"),  # 60%, 50% to a
    ]
    kept = homologs.filter_redundant(QUERY, rows)
    assert [name for name, _ in kept] == ["a", "b"]


def test_filter_redundant_max_hits():
    """At most max_hits hits are kept."""
    rows = [("a", "MKTSYLAKER"), ("b", "LRTAWIGKQE")]
    assert homologs.filter_redundant(QUERY, rows, max_hits=1) == [rows[0]]


def _rows(n):
    """Return n hits that are diverse enough to pass the filter."""
    letters = "ACDEFGHILMNPSTVWY"
    rows = []
    for i in range(n):
        row = list(QUERY)
        for j in range(4):
            row[(i + 3 * j) % len(QUERY)] = letters[(i * 4 + j) % len(letters)]
        rows.append(("hit{0}".format(i), "".join(row)))
    return rows


@pytest.fixture
def no_network(monkeypatch):
    """Make any HTTP request fail the test."""

    def fail(*_args, **_kwargs):
        """Fail on network access."""
        raise AssertionError("unexpected network access")

    monkeypatch.setattr(homologs, "_http", fail)


def test_find_homologs_colabfold(monkeypatch, no_network):
    """Use the ColabFold results when enough hits remain."""
    monkeypatch.setattr(homologs, "colabfold_search", lambda seq: "a3m")
    monkeypatch.setattr(homologs, "parse_a3m", lambda text, seq: _rows(40))
    monkeypatch.setattr(homologs, "filter_redundant", lambda seq, rows: rows[:20])
    kept = homologs.find_homologs(QUERY, email="user@example.org")
    assert len(kept) == 20


def test_find_homologs_falls_back_to_ebi(monkeypatch, no_network):
    """Too few ColabFold hits: the EBI search is used."""
    monkeypatch.setattr(homologs, "colabfold_search", lambda seq: "a3m")
    monkeypatch.setattr(homologs, "parse_a3m", lambda text, seq: _rows(2))
    monkeypatch.setattr(homologs, "ebi_search", lambda seq, email: "xml")
    monkeypatch.setattr(homologs, "parse_ebi_json", lambda text, seq: _rows(30))
    monkeypatch.setattr(homologs, "filter_redundant", lambda seq, rows: rows)
    kept = homologs.find_homologs(QUERY, email="user@example.org")
    assert len(kept) == 30


def test_find_homologs_colabfold_error(monkeypatch, no_network):
    """A failing ColabFold search falls back to the EBI."""

    def broken(seq):
        raise homologs.HomologSearchError("down")

    monkeypatch.setattr(homologs, "colabfold_search", broken)
    monkeypatch.setattr(homologs, "ebi_search", lambda seq, email: "xml")
    monkeypatch.setattr(homologs, "parse_ebi_json", lambda text, seq: _rows(12))
    monkeypatch.setattr(homologs, "filter_redundant", lambda seq, rows: rows)
    assert len(homologs.find_homologs(QUERY, email="user@example.org")) == 12


def test_find_homologs_keeps_the_best_partial_result(monkeypatch, no_network):
    """Without the EBI, a few ColabFold hits are still used."""
    monkeypatch.setattr(homologs, "colabfold_search", lambda seq: "a3m")
    monkeypatch.setattr(homologs, "parse_a3m", lambda text, seq: _rows(5))
    monkeypatch.setattr(homologs, "filter_redundant", lambda seq, rows: rows)
    assert len(homologs.find_homologs(QUERY)) == 5


def test_find_homologs_not_enough(monkeypatch, no_network):
    """Too few homologs from every service: a clear error."""
    monkeypatch.setattr(homologs, "colabfold_search", lambda seq: "a3m")
    monkeypatch.setattr(homologs, "parse_a3m", lambda text, seq: _rows(2))
    monkeypatch.setattr(homologs, "filter_redundant", lambda seq, rows: rows)
    with pytest.raises(homologs.HomologSearchError, match="upload a multiple"):
        homologs.find_homologs(QUERY)


def test_wait_time_limit(monkeypatch):
    """A search that never finishes stops at the time limit."""
    monkeypatch.setattr(homologs, "SEARCH_TIMEOUT", 0)
    monkeypatch.setattr(homologs, "POLL_INTERVAL", 0)
    with pytest.raises(homologs.HomologSearchError, match="did not finish"):
        homologs._wait(lambda: None, "Test")


class FakeServer:
    """Answer _http calls from a list of (URL fragment, response) pairs."""

    def __init__(self, responses):
        """Store the responses, in the order they are requested."""
        self.responses = list(responses)
        self.calls = []

    def __call__(self, url, data=None):
        """Return the next response after checking the URL."""
        self.calls.append((url, data))
        fragment, response = self.responses.pop(0)
        assert fragment in url
        return response


def _a3m_archive(text):
    """Return a ColabFold result archive holding uniref.a3m."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        data = text.encode()
        info = tarfile.TarInfo("uniref.a3m")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


def test_colabfold_search(monkeypatch):
    """Submit, poll until complete and read uniref.a3m from the archive."""
    server = FakeServer(
        [
            ("/ticket/msa", b'{"id": "T1", "status": "PENDING"}'),
            ("/ticket/T1", b'{"id": "T1", "status": "RUNNING"}'),
            ("/ticket/T1", b'{"id": "T1", "status": "COMPLETE"}'),
            ("/result/download/T1", _a3m_archive(A3M)),
        ]
    )
    monkeypatch.setattr(homologs, "_http", server)
    monkeypatch.setattr(homologs, "POLL_INTERVAL", 0)
    assert homologs.colabfold_search(QUERY) == A3M
    assert server.calls[0][1]["mode"] == "all"


def test_colabfold_search_error(monkeypatch):
    """A failed ColabFold job raises a HomologSearchError."""
    server = FakeServer([("/ticket/msa", b'{"id": "T1", "status": "ERROR"}')])
    monkeypatch.setattr(homologs, "_http", server)
    with pytest.raises(homologs.HomologSearchError, match="ERROR"):
        homologs.colabfold_search(QUERY)


def test_ebi_search(monkeypatch):
    """Submit with the contact email, poll and fetch the JSON result."""
    server = FakeServer(
        [
            ("/run", b"JOB1"),
            ("/status/JOB1", b"RUNNING"),
            ("/status/JOB1", b"FINISHED"),
            ("/result/JOB1/json", EBI_JSON.encode()),
        ]
    )
    monkeypatch.setattr(homologs, "_http", server)
    monkeypatch.setattr(homologs, "POLL_INTERVAL", 0)
    assert homologs.ebi_search(QUERY, "user@example.org") == EBI_JSON
    assert server.calls[0][1]["email"] == "user@example.org"


def test_ebi_search_failure(monkeypatch):
    """A failed EBI job raises a HomologSearchError."""
    server = FakeServer([("/run", b"JOB1"), ("/status/JOB1", b"FAILURE")])
    monkeypatch.setattr(homologs, "_http", server)
    with pytest.raises(homologs.HomologSearchError, match="FAILURE"):
        homologs.ebi_search(QUERY, "user@example.org")


class FakeResponse(io.BytesIO):
    """Response object usable as a context manager."""

    def __enter__(self):
        """Return the response itself."""
        return self

    def __exit__(self, *args):
        """Close the response."""
        self.close()


def test_http_retries(monkeypatch):
    """Network errors are retried."""
    attempts = []

    def urlopen(request, timeout):
        attempts.append(request.full_url)
        if len(attempts) < 3:
            raise urllib.error.URLError("unreachable")
        return FakeResponse(b"ok")

    monkeypatch.setattr(homologs.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(homologs, "HTTP_RETRY_WAIT", 0)
    assert homologs._http("https://example.org/x", {"a": 1}) == b"ok"
    assert len(attempts) == 3


def test_http_gives_up(monkeypatch):
    """After the last retry the error is raised as a HomologSearchError."""

    def urlopen(request, timeout):
        raise urllib.error.URLError("unreachable")

    monkeypatch.setattr(homologs.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(homologs, "HTTP_RETRY_WAIT", 0)
    with pytest.raises(homologs.HomologSearchError, match="unreachable"):
        homologs._http("https://example.org/x")


def test_http_only_https():
    """Only https URLs are requested."""
    with pytest.raises(homologs.HomologSearchError, match="https"):
        homologs._http("file:///etc/passwd")
