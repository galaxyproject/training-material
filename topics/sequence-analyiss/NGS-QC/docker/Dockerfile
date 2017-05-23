# Galaxy - NGS QC training material

FROM bgruening/galaxy-ngs-preprocessing

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND NGS QC
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Add a Galaxy tour
# ADD <URL to dedicated tour in https://github.com/galaxyproject/galaxy-tours> $GALAXY_ROOT/config/plugins/tours/ngs-qc.yaml
