# Using WHISCY

After it has been installed according to the instructions in the [INSTALLATION.md](INSTALLATION.md) file, you can use WHISCY to predict the binding site of a protein.

## Example

Below we will give an example of how to use WHISCY to predict the binding site of the protein 1PPE.

### Get an input PDB

Download the PDB file of the protein 1PPE from the [RCSB PDB](https://www.rcsb.org/structure/1PPE) website.

```bash
wget https://files.rcsb.org/download/1PPE.pdb
```

### Setup the prediction

```text
$ whiscy_setup 1PPE.pdb E

  2024-04-22 11:43:39,469 cli_setup:201 INFO - PDB structure with chain E saved to 1PPE_E.pdb
  2024-04-22 11:43:39,494 cli_setup:208 INFO - Atom accessibility calculated to 1PPE_E.rsa
  2024-04-22 11:43:39,495 cli_setup:212 INFO - Surface and buried residues calculated
  2024-04-22 11:43:39,513 cli_setup:255 INFO - HSSP file not found, fallback to generating MSA with blastp
  2024-04-22 11:43:39,513 cli_setup:262 INFO - Running blastp via NCBI
  2024-04-22 11:44:43,926 cli_setup:264 INFO - Result stored in 1PPE_E_blast.xml
  2024-04-22 11:44:43,933 cli_setup:277 INFO - Generating MSA using MUSCLE...
  2024-04-22 11:44:44,027 cli_setup:280 INFO - Done.
  2024-04-22 11:44:44,027 cli_setup:283 INFO - Converting MSA file to Phylseq format...
  2024-04-22 11:44:44,027 cli_setup:292 INFO - 1PPE_E.phylseq file written
  2024-04-22 11:44:44,046 cli_setup:319 INFO - Protdist calculated
  2024-04-22 11:44:44,054 cli_setup:326 INFO - Conversion table file generated
  2024-04-22 11:44:44,054 cli_setup:328 INFO - Whiscy setup finished
```

This script will generate a set of files needed for the prediction step with `whiscy`. Here it is a list of the generated files in our 1PPE complex.

> Follow the links to see a set of pre-generated files for the 1PPE protein.

| File name                                | Explanation                                                        |
| ---------------------------------------- | ------------------------------------------------------------------ |
| [1ppe.hssp](example/1ppe.hssp)           | Multiple sequence alignment download from the HSSP database        |
| [1ppe.hssp.bz2](example/1ppe.hssp.bz2)   | HSSP MSA file compressed                                           |
| [1ppe.pdb](example/1ppe.pdb)             | PDB file download from the Protein Data Bank                       |
| [1ppe_E.pdb](example/1ppe_E.pdb)         | 1ppe.pdb parsed to select only the given `chain_id`                |
| [1ppe_E.rsa](example/1ppe_E.rsa)         | SASA output of `freesasa` in `NACCESS` format of `1ppe_E.pdb` file |
| [1ppe_E.fasta](example/1ppe_E.fasta)     | Sequence of 1ppe_E.pdb. Alternative residues have been removed     |
| [1ppe_E.phylseq](example/1ppe_E.phylseq) | MSA file translated from HSSP to PHYLIP format                     |
| [1ppe_E.conv](example/1ppe_E.conv)       | PDB residue numeration to FASTA sequence numeration                |
| [1ppe_E.out](example/1ppe_E.out)         | Output of the `protdist` software on 1ppe_E.pdb                    |
| [1ppe_E.sur](example/1ppe_E.sur)         | >15 % surface residue list according to `sa_pred_cutoff` cutoff    |
| [1ppe_E.suract](example/1ppe_E.suract)   | >40 % surface residue list according to `sa_act_cutoff` cutoff     |
| [1ppe_E.lac](example/1ppe_E.lac)         | 0-15 % accessible residue list                                     |

### Make the interface prediction

WHISCY needs four input files (generated before):

- `surface_list (.sur)`: List of residues in the interface.
- `conversion_table (.conv)`,represents the mapping of the PDB file residue numeration into the FASTA sequence numeration
- `alignment_file (.phylseq)` is the MSA file in PHYLIP format
- `distance_file (.out)` is the output of Protdist software

```text
$ whiscy 1PPE_E.sur 1PPE_E.conv 1PPE_E.phylseq 1PPE_E.out -o 1PPE_E.cons
  whiscy [INFO] Parsing surface list...
  whiscy [INFO] Loading conversion table...
  whiscy [INFO] Converting...
  whiscy [INFO] Initializing score calculation...
  whiscy [INFO] Calculating scores...
  whiscy [INFO] Subtracting average value ...
  whiscy [INFO] Sorting scores...
  whiscy [INFO] Writing scores...
  whiscy [INFO] Prediction written to 1PPE_E.cons

  "I’m on a whisky diet. I’ve lost three days already."
    -  Tommy Cooper
```

The file `1PPE_E.cons` will contain the initial predictions, but we can improve them by adding the interface propensities and smoothing the surface.

```text
$ head 1PPE_E.cons
  0.10595   N25
  0.10595   N72
  0.10595   N74
  0.10595   N79
  0.10595  N100
  0.10595  N179
  0.10572   N95
  0.10572  N115
  0.10446   S37
  0.10446   S49
```

#### Add Interface propensities

```text
$ whiscy_consadjust 1PPE_E.cons -o 1PPE_E.acons

consadjust [INFO] Reading input files
consadjust [INFO] Subtracting average value...
consadjust [INFO] Sorting scores...
consadjust [INFO] Writing scores...
```

#### Smooth the surface

```text
$ whiscy_resdist 1PPE_E.pdb 1PPE_E.conv 1PPE_E.rd

residue_distance [INFO] Reading conversion table
residue_distance [INFO] Reading PDB structure from 1PPE_E.pdb
residue_distance [INFO] Residue distances written to 1PPE_E.rd
```

```text
$ whiscy_parasmooth 1PPE_E.acons 1PPE_E.cons 1PPE_E.rd -o 1PPE_E.pscons

parasmooth [INFO] Reading input files
parasmooth [INFO] Calculating parameter smoothing
parasmooth [INFO] Result written to 1PPE_E.pscons
```

`1PPE_E.pscons` will contain the final predictions! 🎉

```text
head 1PPE_E.pscons
 0.16475   V75
 0.14702   I73
 0.14633   L114
 0.14547   F82
 0.14297   V76
 0.12907   R117
 0.12127   Q81
 0.11694   N115
 0.11630   N72
 0.11512   N25
```

## Visualizing the results

You can visualize the predictions in the input structure using `whiscy_bfactor`, this script will add the WHISCY scores to the B-factor column of the PDB file.

```text
$ whiscy_bfactor 1PPE_E.pdb 1PPE_E_whiscy.pdb 1PPE_E.pscons

1PPE_E_whiscy.pdb PDB file with WHISCY scores in B-factor column has been created
```

Now open the structure in your favorite molecular viewer and check the WHISCY scores in the B-factor column.
