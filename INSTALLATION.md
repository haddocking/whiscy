# Install Whiscy

Below we will go into a step-by-step guide on how to install Whiscy. Instructions are provived for Ubuntu Linux, in other systems the installation process may vary.

## Third-Party

Whiscy depends on the following third-party software:

- Muscle
- FreeSASA
- HSSPCONV
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

### Troubleshooting

If you get an error such as `ValueError: MUSCLE_BIN not found in system variables`, make sure you have followed the steps in [THIRDPARTY.md](THIRDPARTY.md) to install the third-party dependencies.
