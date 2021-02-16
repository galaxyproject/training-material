
# Galaxy - Using Galaxy tools to generate new Galaxy tools
#
# to build the docker image, go to root of training repo and
#    docker build -t tool-generators -f topics/tool-generators/docker/Dockerfile .
#
# to run image:
#    docker run -p "9090:9090" -t tool-generators
#    use -d to automatically dowload the datalibraries in the container
# ToolFactory planemo will be available on host:9090

FROM ubuntu:latest

MAINTAINER Ross Lazarus
ADD topics/tool-generators/docker/fixplanemo.sh /fixplanemo.sh
ENV TARGDIR "/galaxy-central"
ENV PDIR "/planemo"
RUN chmod a+x /fixplanemo.sh \
&& apt update -y && apt install -y python3-dev gcc build-essential python3-venv python3-wheel nano curl wget git python3-setuptools gnupg curl \
&& curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
&& apt upgrade -y \
&& mkdir -p $TARGDIR \
&& curl -L -s https://github.com/galaxyproject/galaxy/archive/dev.tar.gz | tar xzf - --strip-components=1 -C $TARGDIR \
&& git clone --recursive https://github.com/fubar2/planemo.git $PDIR \
&& cd $PDIR \
&& mkdir mytools \
&& python3 setup.py build \
&& python3 setup.py install \
&& /bin/bash /fixplanemo.sh \
&& apt-get clean && apt-get purge \
&&  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
## fixplanemo adds configparser, planemo does not work in docker otherwise... No clue why but planemo goes all pear shaped.
## pip reports that it's missing - installing it explicitly seems to do some kind of magic
## galaxyxml may not be needed but it was "missing" at one point...

ENV GALAXY_CONFIG_BRAND "GTN: Using Galaxy tools to generate new Galaxy tools"
EXPOSE 9090

ENTRYPOINT ["/usr/local/bin/planemo" ,"tool_factory", "--galaxy_root" ,"/galaxy-central", "--port", "9090", "--host", "0.0.0.0"]

# copy the tutorials directory for your topic
#ADD topics/tool-generators/tutorials/ /tutorials/

# install everything for tutorials
#ADD bin/docker-install-tutorials.sh /setup-tutorials.sh
#ADD bin/mergeyaml.py /mergeyaml.py
#ADD bin/data_libarary_download.sh /data_libarary_download.sh
#RUN cp /galaxy-central/config/datatypes_conf.xml.sample /galaxy-central/config/datatypes_conf.xml \
#&& sed -i 's/<\/registration>/<datatype extension="tgz" type="galaxy.datatypes.binary:Binary" subclass="true" mimetype="multipart\/x-gzip" display_in_upload="true"\/> <\/registration>/' /galaxy-central/config/datatypes_conf.xml \
#&& sed -i 's/<datatype extension="html"/<datatype extension="html" display_in_upload="true"/' /galaxy-central/config/datatypes_conf.xml \
#&& /setup-tutorials.sh

#ENTRYPOINT ["/data_libarary_download.sh"]
