"""Tests for the remote homolog search."""

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

EBI_XML = """<?xml version="1.0"?>
<EBIApplicationResult xmlns="http://www.ebi.ac.uk/schema">
 <SequenceSimilaritySearchResult>
  <hits total="2">
   <hit number="1" database="uniprotkb" id="SP:P1" ac="P1">
    <alignments total="1">
     <alignment number="1">
      <querySeq start="3" end="10">TAY-IAKQR</querySeq>
      <pattern>TAY IAKQR</pattern>
      <matchSeq start="5" end="13">TAYWIAKQK</matchSeq>
     </alignment>
    </alignments>
   </hit>
   <hit number="2" database="uniprotkb" id="TR:Q2" ac="Q2">
    <alignments total="1">
     <alignment number="1">
      <querySeq start="1" end="10">MKTAYIAKQR</querySeq>
      <matchSeq start="1" end="9">MKT-YIGKQR</matchSeq>
     </alignment>
    </alignments>
   </hit>
  </hits>
 </SequenceSimilaritySearchResult>
</EBIApplicationResult>
"""


def test_parse_a3m():
    """Insertions are removed and the query row is skipped."""
    rows = homologs.parse_a3m(A3M, QUERY)
    assert rows == [
        ("UniRef100_A", "MKTAYIAKQK"),
        ("UniRef100_B", "MKTAYIGKQR"),
        ("UniRef100_C", "--TAYLAKQR"),
    ]


def test_parse_ebi_xml():
    """Hits are placed on the query columns; query gaps are dropped."""
    rows = homologs.parse_ebi_xml(EBI_XML, QUERY)
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
    """n hits that are diverse enough to pass the filter."""
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
    """ColabFold results are used when enough hits remain."""
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
    monkeypatch.setattr(homologs, "parse_ebi_xml", lambda text, seq: _rows(30))
    monkeypatch.setattr(homologs, "filter_redundant", lambda seq, rows: rows)
    kept = homologs.find_homologs(QUERY, email="user@example.org")
    assert len(kept) == 30


def test_find_homologs_colabfold_error(monkeypatch, no_network):
    """A failing ColabFold search falls back to the EBI."""

    def broken(seq):
        raise homologs.HomologSearchError("down")

    monkeypatch.setattr(homologs, "colabfold_search", broken)
    monkeypatch.setattr(homologs, "ebi_search", lambda seq, email: "xml")
    monkeypatch.setattr(homologs, "parse_ebi_xml", lambda text, seq: _rows(12))
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
