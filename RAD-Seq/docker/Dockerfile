# Galaxy - RAD-Seq training material

FROM engineson/galaxy-rad-seq

MAINTAINER LE Bras Yvan, yvanlebras & Moreews François, fjrmoreews

ENV GALAXY_CONFIG_BRAND RAD-Seq training
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Add a Galaxy tour
# ADD <URL to dedicated tour in https://github.com/galaxyproject/galaxy-tours> $GALAXY_ROOT/config/plugins/tours/rna-seq.yaml