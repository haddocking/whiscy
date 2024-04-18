#==============================================================================================
FROM python:3.11 AS base

LABEL author="Rodrigo V. Honorato <r.vargashonorato@uu.nl>"

ARG SOFTWARE_PATH=/opt/software

#------------------------------------------------------------------------------------------
# System dependencies
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  build-essential \
  libboost-all-dev \
  && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

#------------------------------------------------------------------------------------------
# Build Muscle
WORKDIR ${SOFTWARE_PATH}/muscle3.8.1551
RUN curl https://drive5.com/muscle/muscle_src_3.8.1551.tar.gz | tar xzv && \
  make
ENV MUSCLE_BIN=${SOFTWARE_PATH}/muscle3.8.1551/muscle

#------------------------------------------------------------------------------------------
# Build hsspconv
WORKDIR ${SOFTWARE_PATH}
RUN wget https://github.com/cmbi/hssp/archive/3.1.5.tar.gz && \
  tar -zxvf 3.1.5.tar.gz && \
  cd hssp-3.1.5 && \
  ./autogen.sh && \
  ./configure && \
  make hsspconv
ENV HSSPCONV_BIN=${SOFTWARE_PATH}/hssp-3.1.5/hsspconv

#------------------------------------------------------------------------------------------
# Build freesasa
WORKDIR ${SOFTWARE_PATH}
RUN  wget https://github.com/mittinatten/freesasa/releases/download/2.0.3/freesasa-2.0.3.tar.gz  && \
  tar -zxvf freesasa-2.0.3.tar.gz && \
  cd freesasa-2.0.3 && \
  ./configure --disable-json --prefix=`pwd` && \
  make && \
  make install
ENV FREESASA_BIN=${SOFTWARE_PATH}/freesasa-2.0.3/bin/freesasa

#------------------------------------------------------------------------------------------
# Install Poetry
RUN pip install --no-cache-dir poetry==1.8.2 \
  && poetry config virtualenvs.create false

# Install Whiscy
WORKDIR ${SOFTWARE_PATH}/whiscy
COPY . .
RUN poetry install --no-dev

#------------------------------------------------------------------------------------------
# Build protdist
WORKDIR ${SOFTWARE_PATH}/whiscy/bin/protdist
RUN sh compile.sh
ENV PROTDIST_BIN=${SOFTWARE_PATH}/whiscy/bin/protdist/protdist

#------------------------------------------------------------------------------------------
# Set data directory

WORKDIR /data

###############################################################################################
# No entrypoint here because Whiscy runs multiple commands
###############################################################################################

#==============================================================================================

FROM base AS test

# RUN pip install pytest coverage pytest pytest-cov hypothesis
WORKDIR ${SOFTWARE_PATH}/whiscy
RUN poetry install


#==============================================================================================

FROM test AS dev

ARG USERNAME=dev
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
  && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
  # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
  && apt-get update \
  && apt-get install -y sudo \
  && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
  && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME