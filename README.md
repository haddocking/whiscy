![WHISCY](media/whiscy_logo.png)

## WHat Information does Surface Conservation Yield?

WHISCY is a program to predict protein-protein interfaces. It is primarily based on conservation, but it also takes into account structural information. A sequence alignment is used to calculate a prediction score for each surface residue of your protein.

This repository contains the Python 3 implementation of the original code developed by [Sjoerd de Vries](https://scholar.google.de/citations?user=fpjNl3wAAAAJ&hl=en) and originally published in [2006](http://dx.doi.org/doi:10.1002/prot.20842):

> *De Vries SJ, van Dijk ADJ, and Bonvin AMJJ*<br>
> WHISCY: What information does surface conservation yield? Application to data-driven docking.<br>
> *Proteins: Struc. Funct. & Bioinformatics*; 2006, **63**(3): 479-489.

## 1. Installation

WHISCY needs the following software to be installed:

* **Python 3** (tested in version 3.6.6)
* **Python 3 libraries**: numpy (1.15.0), nose (1.3.7), biopython (1.71)
* **GNU CC compiler** (gcc, tested on version 6.4.0)
* **FreeSASA** (tested on version 2.0.3)
* **MUSCLE** (only if you want WHISCY to do automotically the multiple sequence alignment for you, tested on version 3.8.31)

Software version is indicative except for Python, which has to be from the 3.6.x series, and freesasa from the series 2.x and above.

### 1.1. Installation in macOS

### 1.1.1. Python3, GCC and libraries
Using [Macports](https://www.macports.org/), you can install Python 3 and the necessary libraries, GCC and git:

```bash
sudo port install gcc6 git python36 py36-biopython py36-numpy py36-nose
```

Similar installation should be possible using [Homebrew](https://brew.sh/) instead of Macports.

If you already have `python3` and `pip3` installed, it is completely OK to install `bio`, `numpy`and `nose` libraries using `pip3`.


### 1.1.2. MUSCLE
To install MUSCLE, go to [the official download site](https://www.drive5.com/muscle/downloads.htm), download the Mac OS X version suitable for your architecture (32 or 64bit) and follow the instructions provided by the authors.

### 1.1.3. FreeSASA
Go to the [FreeSASA main site](https://freesasa.github.io/) and follow the `Quick-start guide`. We don't need the Python bindings as we will be calling freesasa binary from command line.

**Make sure freesasa binary is in your path:**

```bash
$ freesasa --version
FreeSASA 2.0.3
License: MIT <http://opensource.org/licenses/MIT>
If you use this program for research, please cite:
  Simon Mitternacht (2016) FreeSASA: An open source C
  library for solvent accessible surface area calculations.
  F1000Research 5:189.

Report bugs to <https://github.com/mittinatten/freesasa/issues>
Home page: <http://freesasa.github.io>
```

### 1.1.4. WHISCY
Then, the next step is to clone the repository

```bash
git clone https://github.com/haddocking/whiscy.git
cd whiscy
pwd
```

With `pwd`, you will get the directory where you have cloned `whiscy`. Please, copy that directory path because you will need to specify it in your `.bashrc` or `.bash_profile` file:

Edit your `.bashrc` or `.bash_profile` and add the following lines:

```bash
# Whiscy
export WHISCY_PATH=/PATH/TO/WHISCY
export PYTHONPATH=$PYTHONPATH:${WHISCY_PATH}
export WHISCY_BIN=${WHISCY_PATH}/whiscy.py
export PATH=$PATH:${WHISCY_PATH}
```

You have to change `/PATH/TO/WHISCY` according to the directory pointed by the `pwd` command.

Now, we compile the `protdist` software: 

```bash
cd $WHISCY_PATH
cd bin/protdist
./compile.sh
./protdist 
```

If we see an output like this:

```
Too few arguments for this modified version of PROTDIST
Usage: protdist <infile> <outfile>
```
everythin is ready.

There is only one final step where we tell WHISCY where to find MUSCLE binary. Edit `$WHISCY_PATH/etc/local.json`:

```json
{
  "ALIGN": {
    "MUSCLE_BIN": "/path/to/bin/muscle/muscle3.8.31_i86darwin64"
  },
  "CUTOFF": {
    "sa_pred_cutoff": 15.0,
    "sa_act_cutoff": 40.0,
    "air_cutoff": 0.18,
    "air_dist_cutoff": 6.5
  },
  "AIR": {
    "air_pro_percentage": 10.0,
    "air_wm_pro_or": 98.52,
    "air_wm_whis_or": 0.370515,
    "air_wm_pro_and": 55.42,
    "air_wm_whis_and": 0.106667
  }
}
```

Change the `MUSCLE_BIN` variable to the correct path of your MUSCLE binary.



