#==============================================================================================
FROM python:3.6 as base

LABEL author="Rodrigo V. Honorato <r.vargashonorato@uu.nl>"

RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  build-essential \
  libboost-all-dev \
  && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/software

# Install Whiscy
RUN git clone https://github.com/haddocking/whiscy.git && \
  cd whiscy && \
  git checkout tags/v1.1.1 && \
  cd /opt/software/whiscy/bin/protdist && \
  sh compile.sh

RUN pip install biopython==1.79

WORKDIR /opt/software/whiscy/bin

# Build Muscle
WORKDIR /opt/software/whiscy/muscle3.8.1551
RUN curl https://drive5.com/muscle/muscle_src_3.8.1551.tar.gz | tar xzv && \
  make && \
  mv muscle /opt/software/whiscy/bin/muscle3.8.1551 && \
  sed -i "s/\/Users\/bjimenez\/bin\/muscle\/muscle3.8.31_i86darwin64/\/opt\/software\/whiscy\/bin\/muscle3.8.1551/g" /opt/software/whiscy/etc/local.json

# Install hsspconv
RUN wget https://github.com/cmbi/hssp/archive/3.1.5.tar.gz && \
  tar -zxvf 3.1.5.tar.gz && \
  cd hssp-3.1.5 && \
  ./autogen.sh && \
  ./configure && \
  make hsspconv && \
  mv hsspconv ../ && \
  sed -i "s/\/Users\/bjimenez\/bin\/hssp\/hsspconv/\/opt\/software\/whiscy\/bin\/hsspconv/g" /opt/software/whiscy/etc/local.json


# Install freesasa
RUN  wget https://github.com/mittinatten/freesasa/releases/download/2.0.3/freesasa-2.0.3.tar.gz  && \
  tar -zxvf freesasa-2.0.3.tar.gz && \
  cd freesasa-2.0.3 && \
  ./configure --disable-json --prefix=/opt/software/whiscy/bin/freesasa && \
  make && make install

# WHISCY exports
ENV WHISCY_PATH=/opt/software/whiscy
ENV PYTHONPATH="{$PYTHONPATH}:${WHISCY_PATH}"
ENV WHISCY_BIN="${WHISCY_PATH}/whiscy.py"
ENV PATH="${WHISCY_PATH}:${WHISCY_PATH}/bin/freesasa/bin:${PATH}"
#==============================================================================================

FROM base AS test

RUN pip install pytest coverage pytest pytest-cov hypothesis

WORKDIR /opt/software/whiscy

COPY test/ tests/

#==============================================================================================
