# Install Whiscy

Below we will go into a step-by-step guide on how to install Whiscy. Instructions are provived for Ubuntu Linux, in other systems the installation process may vary.

## Third-Party

Whiscy depends on the following third-party software:

- FreeSASA
- Protdist

The installation of these dependencies is described in the [THIRDPARTY.md](THIRDPARTY.md) file.

## Whiscy

It can be easily installed from `pip`

```bash
pip install whiscy
```

Check the installation

```text
$ whiscy -h
usage: whiscy [-h] [-o output_file] [--version]
              surface_list conversion_table alignment_file distance_file

positional arguments:
  surface_list          Surface list
  conversion_table      Conversion table
  alignment_file        Alignment file
  distance_file         Distance file

options:
  -h, --help            show this help message and exit
  -o output_file, --output output_file
                        If set, output prediction to this file
  --version             show program's version number and exit
```

### Environment variables

WHISCY reads these environment variables:

| Variable | Required | Purpose |
| --- | --- | --- |
| `FREESASA_BIN` | yes | Path to the `freesasa` binary |
| `PROTDIST_BIN` | yes | Path to the `protdist` binary |
| `WHISCY_EBI_EMAIL` | no | Contact email for EBI BLAST (see below) |
| `WHISCY_SEARCH_TIMEOUT` | no | Time limit per search service, in s (3600) |

Without an alignment, `whiscy_setup` searches for homologs with ColabFold
MMseqs2, and with EBI BLAST when ColabFold fails or finds fewer than 10
usable homologs. The EBI requires a contact email; without
`WHISCY_EBI_EMAIL` this fallback is skipped. See [USAGE.md](USAGE.md).

```bash
export FREESASA_BIN=$HOME/software/freesasa-2.0.3/bin/freesasa
export PROTDIST_BIN=$HOME/software/protdist/protdist
export WHISCY_EBI_EMAIL=your.name@example.org
```

### Troubleshooting

If you get an error such as
`ValueError: FREESASA_BIN not found in system variables`, make sure you have
followed the steps in [THIRDPARTY.md](THIRDPARTY.md) to install the
third-party dependencies.
