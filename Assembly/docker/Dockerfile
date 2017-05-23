# Galaxy - Assembly training material

FROM bgruening/galaxy-stable

MAINTAINER Simon-Gladman, simon.gladman@unimelb.edu.au

ENV GALAXY_CONFIG_BRAND Assembly

# Install tools (normally not needed, all needed tools must be in the base image)
ADD tools.yml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml

# Data libraries (this is only needed to add data libraries like reference genomes)
ADD library_data.yaml $GALAXY_ROOT/library_data.yaml

# Download training data and populate the data library
RUN startup_lite && \
    sleep 50 && \
    . $GALAXY_VIRTUAL_ENV/bin/activate && \
    python $GALAXY_ROOT/setup_data_libraries.py -i $GALAXY_ROOT/library_data.yaml

# Galaxy tour
#ADD <link to galaxy tour on  https://github.com/galaxyproject/galaxy-tours/> $GALAXY_ROOT/config/plugins/tours/<tour-name>.yaml
