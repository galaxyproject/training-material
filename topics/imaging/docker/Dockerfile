# Galaxy - imaging training material
#
# to build the docker image, go to root of training repo and
#    docker build -t galaxy/imaging-training -f topics/imaging/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t galaxy/imaging-training
#    use -d to automatically dowload the datalibraries in the container

FROM quay.io/bgruening/galaxy:20.05

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Imaging"

# copy the tutorials directory for your topic
ADD topics/imaging/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
RUN /setup-tutorials.sh
