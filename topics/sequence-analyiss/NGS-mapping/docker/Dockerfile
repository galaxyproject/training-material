# Galaxy - NGS mapping training material

FROM bgruening/galaxy-ngs-preprocessing

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND NGS Mapping
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Galaxy tour
ADD <link to galaxy tour on  https://github.com/galaxyproject/galaxy-tours/> \
  $GALAXY_ROOT/config/plugins/tours/<tour-name>.yaml
