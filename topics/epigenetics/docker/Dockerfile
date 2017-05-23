# Galaxy - Epigenetics training material

FROM bgruening/galaxy-stable

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Epigenetics"

## Install everything needed for a tutorial
# Copy the tools
COPY tutorials/methylation-seq/tools.yaml $GALAXY_ROOT/methylation-seq/
# Install the tools
RUN install-tools $GALAXY_ROOT/methylation-seq/tools.yaml
# Copy the workflow
COPY tutorials/methylation-seq/workflow.ga $GALAXY_ROOT/methylation-seq/
# Copy the data libraries
COPY tutorials/methylation-seq/data-libraries.yaml $GALAXY_ROOT/methylation-seq/
# Copy the data manager (if there is one). If there is no data manager, remove the last command in the next RUN
COPY tutorials/methylation-seq/data-manager.yaml $GALAXY_ROOT/methylation-seq/
# Install the workflow, the data libraries and the data manager (if any)
RUN startup_lite && \
    sleep 30 && \
    pip install ephemeris -U && \
    workflow-install --workflow_path $GALAXY_ROOT/methylation-seq/workflow.ga -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \
    setup-data-libraries -i $GALAXY_ROOT/methylation-seq/data-libraries.yaml -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \ 
    run-data-managers -u admin@galaxy.org -p admin --config $GALAXY_ROOT/methylation-seq/data-manager.yaml -g http://localhost:8080
# Copy the tour
COPY tutorials/methylation-seq/tour.yaml $GALAXY_ROOT/config/plugins/tours/methylation-seq.yaml
