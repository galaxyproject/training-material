# Galaxy - ChIP-Seq training material

FROM bgruening/galaxy-stable

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: ChIP-Seq"

## Install everything needed for a tutorial
# Copy the tools
COPY tutorials/tal1-binding-site-identification/tools.yaml $GALAXY_ROOT/tal1-binding-site-identification/
# Install the tools
RUN install-tools $GALAXY_ROOT/tal1-binding-site-identification/tools.yaml
# Copy the workflow
COPY tutorials/tal1-binding-site-identification/workflow.ga $GALAXY_ROOT/tal1-binding-site-identification/
# Copy the data libraries
COPY tutorials/tal1-binding-site-identification/data-libraries.yaml $GALAXY_ROOT/tal1-binding-site-identification/
# Copy the data manager (if there is one). If there is no data manager, remove the last command in the next RUN
COPY tutorials/tal1-binding-site-identification/data-manager.yaml $GALAXY_ROOT/tal1-binding-site-identification/
# Install the workflow, the data libraries and the data manager (if any)
RUN startup_lite && \
    sleep 30 && \
    pip install ephemeris -U && \
    workflow-install --workflow_path $GALAXY_ROOT/tal1-binding-site-identification/workflow.ga -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \
    setup-data-libraries -i $GALAXY_ROOT/tal1-binding-site-identification/data-libraries.yaml -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \ 
    run-data-managers -u admin@galaxy.org -p admin --config $GALAXY_ROOT/tal1-binding-site-identification/data-manager.yaml -g http://localhost:8080
# Copy the tour
COPY tutorials/tal1-binding-site-identification/tour.yaml $GALAXY_ROOT/config/plugins/tours/tal1-binding-site-identification.yaml