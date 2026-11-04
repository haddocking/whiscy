# Using WHISCY

After it has been installed according to the instructions in the [INSTALLATION.md](INSTALLATION.md) file, you can use WHISCY to predict the binding site of a protein.

## Example

Below we will give an example of how to use WHISCY to predict the binding site
of trypsin, chain E of the trypsin-inhibitor complex 1PPE. The `example` directory holds
the files of two runs, made with the commands below:

- `example/colabfold`: the homologs are searched for automatically (the
  default); see its [report](example/colabfold/report.html) and
  [log](example/colabfold/whiscy.log);
- `example/hssp`: the alignment comes from an HSSP file; see its
  [report](example/hssp/report.html) and [log](example/hssp/whiscy.log).

### Get an input PDB

Download the PDB file of the protein 1PPE from the [RCSB PDB](https://www.rcsb.org/structure/1PPE) website.

```bash
wget https://files.rcsb.org/download/1PPE.pdb -O 1ppe.pdb
```

### Setup the prediction

`whiscy_setup` extracts the chain, calculates the accessibility of its
residues, and makes the alignment WHISCY needs. The alignment comes from one
of:

- an alignment file of your own, given with `--alignment` (and
  `--alignment_format`, by default `fasta`);
- an HSSP file, given with `--hssp <name>` for `<name>.hssp` (see below);
- a search for homologous sequences, when neither is given.

#### Searching for homologs (example/colabfold)

```text
$ whiscy_setup 1ppe.pdb E

cli_setup:100 INFO - WHISCY 1.6.0: whiscy_setup 1ppe.pdb E
cli_setup:153 INFO - PDB structure with chain E saved to 1ppe_E.pdb
cli_setup:160 INFO - Atom accessibility calculated to 1ppe_E.rsa
cli_setup:164 INFO - Surface and buried residues calculated
cli_setup:186 INFO - HSSP file not found, searching for homologous sequences
homologs:233 WARNING - No contact email set, the EBI BLAST fallback is disabled
homologs:237 INFO - Searching homologs with ColabFold MMseqs2
homologs:251 INFO - ColabFold MMseqs2: 5000 hits, 250 kept (redundant ones removed)
cli_setup:205 INFO - Homologs stored in 1ppe_E_homologs.fasta
cli_setup:212 INFO - 1ppe_E.phylseq written with 250 homologs
cli_setup:251 INFO - Alignment mapped: 220 of 220 PDB residues, 220 alignment columns
cli_setup:260 INFO - Protdist calculated
cli_setup:267 INFO - Conversion table file generated
cli_setup:269 INFO - Whiscy setup finished
```

The log is shortened here: each line also starts with the date and time.

Without an alignment (`--alignment`), `whiscy_setup` searches for homologous
sequences of the structure's sequence:

1. With the [ColabFold MMseqs2 server](https://github.com/sokrypton/ColabFold)
   (UniRef). This usually takes seconds to a few minutes.
2. If that fails or finds fewer than 10 usable homologs, with
   [EBI BLAST](https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast) against
   UniProtKB. The EBI requires a contact email, given with the
   `WHISCY_EBI_EMAIL` environment variable; without it this step is skipped.

Near-identical hits carry no conservation signal, so the hits are filtered:
at most 90% identity to the query and to each other, at least 50% of the
query covered, up to 250 hits. The result is stored in
`<name>_<chain>_homologs.fasta` and reused if the command is run again. If
fewer than 3 homologs are found, `whiscy_setup` stops and asks for an
alignment. `WHISCY_SEARCH_TIMEOUT` sets the time limit per service in seconds
(default 3600). The search results change as the sequence databases grow, so
a new search can give other homologs than the ones in
[1ppe_E_homologs.fasta](example/colabfold/1ppe_E_homologs.fasta).

#### From an HSSP file (example/hssp)

```text
$ whiscy_setup 1ppe.pdb E --hssp 1ppe

cli_setup:100 INFO - WHISCY 1.6.0: whiscy_setup 1ppe.pdb E --hssp 1ppe
cli_setup:153 INFO - PDB structure with chain E saved to 1ppe_E.pdb
cli_setup:160 INFO - Atom accessibility calculated to 1ppe_E.rsa
cli_setup:164 INFO - Surface and buried residues calculated
cli_setup:177 INFO - HSSP file converted to PHYLIP format
cli_setup:251 INFO - Alignment mapped: 220 of 220 PDB residues, 223 alignment columns
cli_setup:260 INFO - Protdist calculated
cli_setup:267 INFO - Conversion table file generated
cli_setup:269 INFO - Whiscy setup finished
```

HSSP files can be obtained from the
[XSSP server](https://www3.cmbi.umcn.nl/xssp/) of the CMBI: choose the output
"HSSP" and give a PDB code, or use its
[REST API](https://www3.cmbi.umcn.nl/xssp/api/examples), e.g.:

```bash
# Submit the job, then poll its status and fetch the result
curl -X POST -F "data=1ppe" https://www3.cmbi.umcn.nl/xssp/api/create/pdb_id/hssp_hssp/
curl https://www3.cmbi.umcn.nl/xssp/api/status/pdb_id/hssp_hssp/<id>/
curl https://www3.cmbi.umcn.nl/xssp/api/result/pdb_id/hssp_hssp/<id>/
```

The result is JSON; its `result` field is the HSSP file. For 1PPE this is the
same alignment as [1ppe.hssp](example/hssp/1ppe.hssp) (generated in 2017, 330
sequences). WHISCY reads this classic HSSP format, not the Stockholm one.

#### Files

`whiscy_setup` writes the files needed for the prediction step with
`whiscy`:

| File name | Explanation |
| --- | --- |
| [1ppe_E.pdb](example/colabfold/1ppe_E.pdb) | 1ppe.pdb parsed to select only the given `chain_id` |
| [1ppe_E.rsa](example/colabfold/1ppe_E.rsa) | SASA output of `freesasa` in `NACCESS` format of `1ppe_E.pdb` file |
| [1ppe_E.fasta](example/colabfold/1ppe_E.fasta) | Sequence of 1ppe_E.pdb. Alternative residues have been removed |
| [1ppe_E_homologs.fasta](example/colabfold/1ppe_E_homologs.fasta) | Homologs found by the search (not with `--hssp` or `--alignment`) |
| [1ppe_E.phylseq](example/colabfold/1ppe_E.phylseq) | Alignment in PHYLIP format |
| [1ppe_E.conv](example/colabfold/1ppe_E.conv) | PDB residue numeration to alignment numeration |
| [1ppe_E.out](example/colabfold/1ppe_E.out) | Output of the `protdist` software on the alignment |
| [1ppe_E.sur](example/colabfold/1ppe_E.sur) | >15 % surface residue list according to `sa_pred_cutoff` cutoff |
| [1ppe_E.suract](example/colabfold/1ppe_E.suract) | >40 % surface residue list according to `sa_act_cutoff` cutoff |
| [1ppe_E.lac](example/colabfold/1ppe_E.lac) | 0-15 % accessible residue list |

### Make the interface prediction

WHISCY needs four input files (generated before):

- `surface_list (.sur)`: List of residues in the interface.
- `conversion_table (.conv)`,represents the mapping of the PDB file residue numeration into the FASTA sequence numeration
- `alignment_file (.phylseq)` is the MSA file in PHYLIP format
- `distance_file (.out)` is the output of Protdist software

```text
$ whiscy 1ppe_E.sur 1ppe_E.conv 1ppe_E.phylseq 1ppe_E.out -o 1ppe_E.cons
  whiscy [INFO] WHISCY 1.6.0: whiscy 1ppe_E.sur 1ppe_E.conv 1ppe_E.phylseq 1ppe_E.out -o 1ppe_E.cons
  whiscy [INFO] Parsing surface list...
  whiscy [INFO] Loading conversion table...
  whiscy [INFO] Converting...
  whiscy [INFO] Initializing score calculation...
  whiscy [INFO] Calculating scores...
  whiscy [INFO] Subtracting average value ...
  whiscy [INFO] Sorting scores...
  whiscy [INFO] Writing scores...
  whiscy [INFO] Prediction written to 1ppe_E.cons

  "I’m on a whisky diet. I’ve lost three days already."
    -  Tommy Cooper
```

The file `1ppe_E.cons` will contain the initial predictions, but we can improve them by adding the interface propensities and smoothing the surface.

```text
$ head -5 1ppe_E.cons
0.28571  N245
0.28489  S214
0.28370  S195
0.26034  Q192
0.25988  I176
```

The same command on the residues of low accessibility (`.lac`) gives their
scores, used for the smoothing below:

```text
whiscy 1ppe_E.lac 1ppe_E.conv 1ppe_E.phylseq 1ppe_E.out -o 1ppe_E.lcons
```

#### Add Interface propensities

```text
$ whiscy_consadjust 1ppe_E.cons -o 1ppe_E.acons

consadjust [INFO] WHISCY 1.6.0: whiscy_consadjust 1ppe_E.cons -o 1ppe_E.acons
consadjust [INFO] Reading input files
consadjust [INFO] Subtracting average value...
consadjust [INFO] Sorting scores...
consadjust [INFO] Writing scores...
```

#### Smooth the surface

```text
$ whiscy_resdist 1ppe_E.pdb 1ppe_E.conv 1ppe_E.rd

residue_distance [INFO] WHISCY 1.6.0: whiscy_resdist 1ppe_E.pdb 1ppe_E.conv 1ppe_E.rd
residue_distance [INFO] Reading conversion table
residue_distance [INFO] Reading PDB structure from 1ppe_E.pdb
residue_distance [INFO] Residue distances written to 1ppe_E.rd
```

```text
$ whiscy_parasmooth 1ppe_E.acons 1ppe_E.lcons 1ppe_E.rd -o 1ppe_E.pscons

parasmooth [INFO] WHISCY 1.6.0: whiscy_parasmooth 1ppe_E.acons 1ppe_E.lcons 1ppe_E.rd -o 1ppe_E.pscons
parasmooth [INFO] Reading input files
parasmooth [INFO] Calculating parameter smoothing
parasmooth [INFO] Result written to 1ppe_E.pscons
```

`1ppe_E.pscons` will contain the final predictions! 🎉 Here for the two
examples (left: homolog search, right: HSSP file):

```text
$ head 1ppe_E.pscons
 0.19753   Y39         0.31595   Y39
 0.17385   F41         0.31046   H57
 0.16250   I73         0.24401   C191
 0.12739   C191        0.23303   Y94
 0.12673   H40         0.23289   I73
 0.12072   Y94         0.21471   F41
 0.10342   H57         0.20514   S195
 0.09771   V75         0.19676   G193
 0.09197   Q192        0.18962   Q192
 0.08481   F82         0.17610   S214
```

Both put the inhibitor-binding site of trypsin at the top (Y39-F41, the
catalytic H57 and S195, C191-G193), but the scores differ with the alignment.

## Visualizing the results

You can visualize the predictions in the input structure using `whiscy_bfactor`, this script will add the WHISCY scores to the B-factor column of the PDB file.

```text
$ whiscy_bfactor 1ppe_E.pdb 1ppe_E_whiscy.pdb 1ppe_E.pscons

whiscy_bfactor [INFO] WHISCY 1.6.0: whiscy_bfactor 1ppe_E.pdb 1ppe_E_whiscy.pdb 1ppe_E.pscons
whiscy_bfactor [INFO] 1ppe_E_whiscy.pdb PDB file with WHISCY scores in B-factor column has been created
```

Now open the structure in your favorite molecular viewer and check the WHISCY scores in the B-factor column.

## Residues for HADDOCK and the HTML report

`whiscy_haddock` selects active and passive residues for
[HADDOCK](https://wenmr.science.uu.nl/haddock2.4/) in two ways:

- active residues with a WHISCY score of at least 0.18: `active.txt`,
  `passive.txt` and `haddock-cutoff.pdb`;
- the 10 best-scored residues as active: `active_top10.txt`,
  `passive_top10.txt` and `haddock-top10.pdb`.

The passive residues are the surface residues next to an active one. The
`.txt` files list the residue numbers, comma-separated; the PDB files have
B-factors of 100 for active, -100 for passive and 0 for the other residues.

```text
$ whiscy_haddock 1ppe_E.pdb 1ppe_E.pscons --surface_list 1ppe_E.sur

whiscy_haddock [INFO] WHISCY 1.6.0: whiscy_haddock 1ppe_E.pdb 1ppe_E.pscons --surface_list 1ppe_E.sur
whiscy_haddock [INFO] Report written to /path/to/report.html
```

In the examples, the active residues are:

| Example | Score ≥ 0.18 | Top 10 |
| --- | --- | --- |
| [Homolog search](example/colabfold/report.html) | COLABFOLD_CUTOFF | COLABFOLD_TOP10 |
| [HSSP file](example/hssp/report.html) | 39, 41, 57, 73, 94, 191, 192, 193, 195 | the same and 214 |

It also writes `report.html` and `state.mvsj` (the MolViewSpec state of the
Mol* view) in the current directory; use `--no-report` to turn this off. The
report shows the structure colored by WHISCY score and by the active and
passive residues of each selection, the residue lists, and the table of
scores. It loads the files of the directory, so open it through a web server,
e.g. `python -m http.server`, then <http://localhost:8000/report.html>.
