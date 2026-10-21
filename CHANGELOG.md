# Changelog

All notable changes to WHISCY, one line per pull request, newest first.
Releases are tagged in git and published as Docker images at
`ghcr.io/haddocking/whiscy:<tag>`.

## [Unreleased]

- Replace the NCBI BLAST search with the ColabFold MMseqs2 server, falling back
  to EBI BLAST (UniProtKB), and filter out redundant hits; drop MUSCLE (#40,
  closes #39)

## [v1.3.1] - 2026-10-20

- Take the HSSP `MASTER` row from the HSSP alignment, so residues with
  insertion codes no longer shift or break the alignment (#38, closes #37)

## [v1.3.0] - 2026-10-20

- Map alignment columns onto the target residues in `whiscy_setup`, so
  conservation scores are assigned to the right residues for uploaded and BLAST
  alignments; pin the Docker base image to `python:3.11-bookworm` (protdist no
  longer compiles on Debian trixie); run the setup regression tests offline;
  exclude tests from Codacy's Bandit check (#35, closes #34)

## [v1.2.1] - 2024-06-27

- Exit with code 1 instead of 0 when a step fails (#33, closes #29, #30, #31,
  #32)
- Integrate with Codacy (#28)
- Update the installation instructions (#26)

## [v1.2.0] - 2024-04-22

- Repository-wide cleanup and linting (#25)
- Update the installation and usage documentation and the README (#24)
- Reorganize the repository (#23, closes #22)
- Add `pyproject.toml` (#21, closes #20)
- Configure third-party binaries through environment variables instead of
  `local.json` (#19, closes #18)
- Update the base Python version to 3.11 (#17, closes #14)

## [v1.1.3] - 2024-04-17

- Fall back to generating the alignment with BLAST (NCBI) and MUSCLE when no
  HSSP file is available (#16)

## [v1.1.2] - 2024-04-15

- Refactor the tests (#15)
- Add a devcontainer configuration (#12)
- Add CI/CD (#11)

## [v1.1.1] - 2024-03-14

- Fix an unbound variable in `get_from_ftp` (#10, closes #9)

## [v1.1.0] - 2021-04-06

- Skip non-standard residues instead of failing (#7, closes #6)
- Support HSSP3 files, converted with `hsspconv` (#5)
- Add the license (#4)
- Link DOIs to the preferred resolver (#3)

## [1.0] - 2018-09-03

- Support older Python versions, warn when configured paths are missing, and
  make protdist compile with older GCC (#2)
