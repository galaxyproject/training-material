
# Galaxy - Translatomics
#
# to build the docker image, go to root of training repo and
#    docker build -t translatomics -f topics/translatomics/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t translatomics
#    use -d to automatically dowload the datalibraries in the container

FROM bgruening/galaxy-stable:latest

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Translatomics"

# copy the tutorials directory for your topic
ADD topics/translatomics/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
ADD bin/data_libarary_download.sh /data_libarary_download.sh
RUN /setup-tutorials.sh

ENTRYPOINT ["/data_libarary_download.sh"]