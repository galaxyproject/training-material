# Galaxy - Introduction training material

FROM bgruening/galaxy-stable

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND Galaxy Introduction
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Add a Galaxy tour
ADD https://raw.githubusercontent.com/galaxyproject/galaxy-tours/master/advanced-galaxy-introduction/Introduction.yaml $GALAXY_ROOT/config/plugins/tours/advanced-galaxy-introduction.yaml
