# Galaxy - DNA Methylation data analysis training material

FROM bgruening/galaxy-epigenetics

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND DNA Methylation
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Install tools
ADD methylc-seq.yaml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml

# Add a Galaxy tour
# ADD <URL to dedicated tour in https://github.com/galaxyproject/galaxy-tours> $GALAXY_ROOT/config/plugins/tours/methylc-seq.yaml
