# Galaxy - <topic name> training material

FROM <base-galaxy-flavor>

MAINTAINER <maintainer-name>, <maintainer-email>

ENV GALAXY_CONFIG_BRAND Topic Name
ENV ENABLE_TTS_INSTALL True

# Enable Conda dependency resolution
ENV GALAXY_CONFIG_CONDA_AUTO_INSTALL=True \
    GALAXY_CONFIG_CONDA_AUTO_INIT=True

# Install tools (normally not needed, all needed tools must be in the base image)
ADD tools.yml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml

# Data libraries (this is only needed to add data libraries like reference genomes)
ADD setup_data_libraries.py $GALAXY_ROOT/setup_data_libraries.py
ADD library_data.yaml $GALAXY_ROOT/library_data.yaml

# Download training data and populate the data library
RUN startup_lite && \
    sleep 50 && \
    . $GALAXY_VIRTUAL_ENV/bin/activate && \
    python $GALAXY_ROOT/setup_data_libraries.py -i $GALAXY_ROOT/library_data.yaml

# Galaxy tour
ADD <link to galaxy tour on  https://github.com/galaxyproject/galaxy-tours/> \
  $GALAXY_ROOT/config/plugins/tours/<tour-name>.yaml
