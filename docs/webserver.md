# Web server manual

![WHISCY](media/whiscy_logo.png)

The WHISCY web server at <https://wenmr.science.uu.nl/whiscy> predicts the
protein-protein interface of one protein chain and gives the active and passive
residues to use as restraints in [HADDOCK](https://www.bonvinlab.org/software/).
You need an account on the portal to submit a job.

## How WHISCY works

WHISCY is primarily based on conservation, but it also takes into account
structural information:

1. WHISCY finds the master sequence: the sequence of the alignment that best
   matches the structure.
2. It estimates the sequence distance (amount of mutation) between the master
   sequence and every other sequence. This sets the amount of mutation expected
   at each position.
3. For each surface residue, it compares the expected with the observed
   mutation. Less change than expected means conservation, and a positive
   WHISCY score.
4. **Interface propensities** adjust the scores: phenylalanines, for example,
   are often in protein-protein interfaces and get a higher score; lysines
   rarely are and get a lower score.
5. **Surface smoothing** spreads the scores over the surface of the structure:
   interfaces form patches, so the neighbours of an interface residue are often
   interface residues too.

## Input

### Structure

The structure of your protein as a PDB file (`.pdb`, at most 10 MB).

### Chain

The chain to predict on. Many PDB files contain several molecules (chains);
WHISCY predicts on one. The chain must be in the `ATOM` records of the file.

### Alignment (optional, recommended)

A multiple sequence alignment of your protein with homologous sequences, in one
of these formats, recognised by the file extension:

| Extension | Format |
| --- | --- |
| `.fasta` | FASTA alignment (default output of MUSCLE) |
| `.aln` | CLUSTAL |
| `.msf` | MSF |
| `.phylseq` | PHYLIP |

- The **first sequence** of the alignment must be the sequence of your protein.
  This can be the sequence of the PDB structure, but only if it is not
  engineered. The structure may cover only part of that sequence.
- A FASTA file must be aligned: at least 2 sequences, all of the same length
  (with gaps).
- WHISCY is robust to a bad alignment quality, so manual correction of the
  alignment is not essential.

Without an alignment, WHISCY builds one from the homologs found by the
[ColabFold MMseqs2](https://github.com/sokrypton/ColabFold) server or, if it is
unavailable, by [EBI BLAST](https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast) against
UniProtKB. A job with an alignment takes a few seconds; without one, it usually
takes a few minutes.

### Options

- **Use interface propensities** (on by default).
- **Use surface smoothing** (on by default). Smoothing requires interface
  propensities.

## Output

The results are shown online and can be downloaded for **14 days**. WHISCY
gives a prediction score to every surface residue that can be matched to the
alignment: the higher the score, the more likely the residue is in a
protein-protein interface.

### Report

The results page shows the WHISCY report (`report.html`):

- a 3D viewer with three views: the surface coloured by WHISCY score, from
  blue (lowest) to white (0) to red (highest), and, with **Next**, the HADDOCK
  residues of each selection (score cutoff, then top 10), active in red and
  passive in green;
- the active and passive residues for HADDOCK, for both selections;
- a table of the scores of all residues, which can be sorted by clicking a
  column header.

For an example, see the [example report](example/report.html) of the 1PPE
chain E prediction.

### Residues for HADDOCK

WHISCY gives two selections of active residues:

- **Score cutoff** (`active.txt`): residues with a WHISCY score of at least
  0.18. Some proteins get few or no active residues this way.
- **Top 10** (`active_top10.txt`): the 10 best-scored residues.

For each, the **passive residues** (`passive.txt`, `passive_top10.txt`) are
the other surface residues with an atom within 6.5 Å of the Cα of an active
residue. The files are comma-separated lists of residue numbers, ready to
paste in HADDOCK.

### Files

The download (`whiscy.zip`) contains all the files of the job, among which:

| File | Content |
| --- | --- |
| `report.html`, `state.mvsj` | The report and its 3D views |
| `<name>.pscons` | The WHISCY scores, one residue per line (see below) |
| `<name>_bfactor.pdb` | The structure, WHISCY scores as B-factors |
| `haddock-cutoff.pdb`, `haddock-top10.pdb` | HADDOCK residues as B-factors |
| `active.txt`, `passive.txt` | The residues for HADDOCK (score cutoff) |
| `active_top10.txt`, `passive_top10.txt` | The residues for HADDOCK (top 10) |
| `whiscy.log` | The log of the job |

Without surface smoothing the scores file is `<name>.acons`, and without
either option `<name>.cons`.

To view `report.html` from the download, serve its directory with
`python -m http.server -d <directory>` and open
<http://localhost:8000/report.html>: browsers do not let a page opened from
disk load the files of the 3D viewer.

### If a job fails

The results page shows the error output of the job, which often says why it
failed (for example an alignment whose first sequence does not match the
structure). For help, ask at <https://ask.bioexcel.eu>.

## Interpreting the predictions

You should not trust WHISCY predictions blindly, but always try to relate them
to information that you already have.

- WHISCY works very well for proteins with a single conserved protein binding
  site.
- It should not be used for non-conserved interactions, such as
  antibody-antigen.
- If you know that your protein has multiple binding sites, verify whether a
  predicted site could be the one you are interested in.

## Running WHISCY locally

WHISCY is also available as a stand-alone program: see
[Installation](INSTALLATION.md) and [Usage](USAGE.md).

## References

To cite WHISCY, please refer to:

- S.J. de Vries, A.D.J. van Dijk and A.M.J.J. Bonvin.
  [WHISCY: What information does surface conservation yield? Application to data-driven docking.](https://doi.org/10.1002/prot.20842)
  *Proteins: Struc. Funct. & Bioinformatics*, **63**, 479-489 (2006).
- R.V. Honorato, P.I. Koukos, B. Jimenez-Garcia, A. Tsaregorodtsev, M. Verlato,
  A. Giachetti, A. Rosato and A.M.J.J. Bonvin.
  [Structural biology in the clouds: The WeNMR-EOSC Ecosystem.](https://doi.org/10.3389/fmolb.2021.729513)
  *Frontiers Mol. Biosci.*, **8**, 729513 (2021).
