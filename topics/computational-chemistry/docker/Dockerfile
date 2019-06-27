
# Galaxy - Computational chemistry
#
# to build the docker image, go to root of training repo and
#    docker build -t computational-chemistry -f topics/computational-chemistry/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t computational-chemistry

FROM bgruening/galaxy-stable

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Computational chemistry"

# prerequisites
RUN pip install ephemeris -U
ADD bin/galaxy-sleep.py /galaxy-sleep.py

# copy the tutorials directory for your topic
ADD topics/computational-chemistry/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
RUN /setup-tutorials.sh