# Galaxy - <topic name> training material

FROM bgruening/galaxy-stable

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Topic Name"

## Install everything needed for a tutorial
# Copy the tools
COPY tutorials/tutorial1/tools.yaml $GALAXY_ROOT/tutorial1/
# Install the tools
RUN install-tools $GALAXY_ROOT/tutorial1/tools.yaml
# Copy the workflow
COPY tutorials/tutorial1/workflow.ga $GALAXY_ROOT/tutorial1/
# Copy the data libraries
COPY tutorials/tutorial1/data-libraries.yaml $GALAXY_ROOT/tutorial1/
# Copy the data manager (if there is one). If there is no data manager, remove the last command in the next RUN
COPY tutorials/tutorial1/data-manager.yaml $GALAXY_ROOT/tutorial1/
# Install the workflow, the data libraries and the data manager (if any)
RUN startup_lite && \
    sleep 30 && \
    pip install ephemeris -U && \
    workflow-install --workflow_path $GALAXY_ROOT/tutorial1/workflow.ga -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \
    setup-data-libraries -i $GALAXY_ROOT/tutorial1/data-libraries.yaml -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD && \ 
    run-data-managers -u admin@galaxy.org -p admin --config $GALAXY_ROOT/tutorial1/data-manager.yaml -g http://localhost:8080
# Copy the tour
COPY tutorials/tutorial1/tour.yaml $GALAXY_ROOT/config/plugins/tours/tutorial1.yaml
