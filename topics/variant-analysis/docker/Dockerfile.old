# Galaxy - Exome sequencing training material

FROM bgruening/galaxy-exome-seq

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND Exome sequencing
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Data libraries
ADD setup_data_libraries.py $GALAXY_ROOT/setup_data_libraries.py
ADD library_data.yaml $GALAXY_ROOT/library_data.yaml

# Download training data and populate the data library
RUN startup_lite && \
    sleep 50 && \
    . $GALAXY_VIRTUAL_ENV/bin/activate && \
    python $GALAXY_ROOT/setup_data_libraries.py -i $GALAXY_ROOT/library_data.yaml

ADD https://raw.githubusercontent.com/galaxyproject/galaxy-tours/8bbd07bd906bf8182f1aa3bc6a573a98bd888667/exom-seq/exom-seq.yaml \
  $GALAXY_ROOT/config/plugins/tours/exom-seq.yaml
