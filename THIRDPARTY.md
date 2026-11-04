# Installation of Third-Party dependencies

These instructions assuming you are installing locally in a Ubuntu Linux system, different steps may be required for other systems.

All will be installed in a `software` directory in your `$HOME` directory.

## System dependendencies

```bash
$ sudo apt-get update && \
  sudo apt-get install -y build-essential
```

## Freesasa

```bash
mkdir -p $HOME/software && cd $HOME/software
wget https://github.com/mittinatten/freesasa/releases/download/2.0.3/freesasa-2.0.3.tar.gz
tar -zxf freesasa-2.0.3.tar.gz && rm freesasa-2.0.3.tar.gz
cd freesasa-2.0.3
./configure --disable-json --prefix=`pwd`
make
make install
export FREESASA_BIN=$HOME/software/freesasa-2.0.3/bin/freesasa
```

## Protdist

Protdist is distributed together with WHISCY, you can find it in the `whiscy` directory. We are working on a better way to install this dependency 🙂

```bash
mkdir -p $HOME/software && cd $HOME/software
git clone https://github.com/haddocking/whiscy
mv whiscy/bin/protdist . && rm -rf whiscy
cd protdist
bash compile.sh
export PROTDIST_BIN=$HOME/software/protdist/protdist
```

---

In the end the system variables that define the third-party dependencies should look like this:

```bash
export FREESASA_BIN=$HOME/software/freesasa-2.0.3/bin/freesasa
export PROTDIST_BIN=$HOME/software/protdist/protdist
```
