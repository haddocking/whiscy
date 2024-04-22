# WHISCY - WHat Information does Surface Conservation Yield?

[![ci](https://github.com/haddocking/whiscy/actions/workflows/test.yml/badge.svg)](https://github.com/haddocking/whiscy/actions/workflows/test.yml)

![WHISCY](media/whiscy_logo.png)

**WHISCY is a program to predict protein-protein interfaces.**

It is primarily based on conservation, but it also takes into account structural information. A sequence alignment is used to calculate a prediction score for each surface residue of your protein.

This repository contains an updated version of WHISCY. The original code was published in the following paper:

- de Vries, S. J., van Dijk, A. D. J. & Bonvin, A. M. J. J. [WHISCY: What information does surface conservation yield? Application to data‐driven docking.](https://doi.org/doi:10.1002/prot.20842) Proteins: Structure, Function, and Bioinformatics vol. 63 479–489 (2006).


## Table of contents

- [How does WHISCY work?](#how-does-whiscy-work)
- [Installation](#installation)
- [Usage](#usage)

## How does WHISCY work?

WHISCY requires a protein structure and a sequence alignment. First, it identifies a master sequence, the sequence that best matches the structure.

The sequence distance (amount of mutation) between the master sequence and all sequences is estimated. This determines the amount of expected mutation.

Then, for each residue, the expected mutation is compared with the observed mutation. Less change than expected means conservation, translated into a positive WHISCY score.

![compare](media/compare.jpg)

Next, the interface propensity is taken into account.

Phenylalanines, for example, are likely to be in a protein-protein interface, so all phenylalanines receive a higher score. Lysines are much less likely to be in a protein-protein interface, so lysines receive a lower score.

Finally, all scores are smoothed over the surface of the protein structure.

Interfaces often form patches, so that neighbours of interface residues often are interface residues, too. The smoothing means that the scores of these neighbours are taken into account.

![interface](media/interface.jpg)

## Installation

Please refer to [INSTALL.md](INSTALL.md) for a detailed guide on how to install WHISCY in macOS and GNU/Linux.

## Usage

Check [USAGE.md](USAGE.md) to learn how to use WHISCY to predict protein-protein interfaces.

