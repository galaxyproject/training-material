
# Galaxy - Materials
#
# to build the docker image, go to root of training repo and
#    docker build -t materials-science -f topics/materials-science/docker/Dockerfile .
#
# to run image:
#    docker run -p "8080:80" -t materials-science
#    use -d to automatically dowload the datalibraries in the container

FROM quay.io/bgruening/galaxy:20.05

MAINTAINER Galaxy Training Material

ENV GALAXY_CONFIG_BRAND "GTN: Materials"

# copy the tutorials directory for your topic
ADD topics/materials-science/tutorials/ /tutorials/

# install everything for tutorials
ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
ADD bin/mergeyaml.py /mergeyaml.py
ADD bin/data_libarary_download.sh /data_libarary_download.sh
RUN /setup-tutorials.sh

ENTRYPOINT ["/data_libarary_download.sh"]