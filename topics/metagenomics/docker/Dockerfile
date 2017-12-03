# Galaxy - metagenomics training material
#
# to build the docker image, go to root of training repo and
#    docker build -t <your_tag> -f topics/<your_topic>/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t <your_tag>

FROM bgruening/galaxy-stable

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Metagenomics"

# copy the tutorials directory for your topic
ADD topics/metagenomics/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
RUN /setup-tutorials.sh
