
# Galaxy - Using Galaxy tools to generate new Galaxy tools
#
# to build the docker image, go to root of training repo and
#    docker build -t tool-generators -f topics/tool-generators/docker/Dockerfile .
#
# to run image:
#    docker run -p "9090:9090" -t tool-generators
# ToolFactory planemo will be available on localhost:9090

FROM ubuntu:latest

MAINTAINER Ross Lazarus
ENV TARGDIR "/galaxy-central"
ENV PDIR "/planemo"
RUN apt update -y && apt install -y -qq python3-dev gcc build-essential python3-venv python3-wheel nano curl wget git python3-setuptools gnupg curl mercurial \
&& curl -sS https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
&& apt upgrade -y \
&& mkdir -p $TARGDIR \
&& curl -L -s https://github.com/galaxyproject/galaxy/archive/dev.tar.gz | tar xzf - --strip-components=1 -C $TARGDIR \
&& git clone --recursive https://github.com/fubar2/planemo.git $PDIR \
&& cd $PDIR \
&& mkdir mytools \
&& python3 setup.py build \
&& python3 setup.py install \
&& planemo conda_init --conda_prefix $PDIR/con \
&& hg clone https://fubar@toolshed.g2.bx.psu.edu/repos/fubar/tacrev  /planemo/tacrev \
&& planemo test --galaxy_root $TARGDIR /planemo/tacrev \
&& cp /galaxy-central/config/datatypes_conf.xml.sample /galaxy-central/config/datatypes_conf.xml \
&& sed -i 's/<\/registration>/<datatype extension="tgz" type="galaxy.datatypes.binary:Binary" subclass="true" mimetype="multipart\/x-gzip" display_in_upload="true"\/> <\/registration>/' /galaxy-central/config/datatypes_conf.xml \
&& sed -i 's/<datatype extension="html"/<datatype extension="html" display_in_upload="true"/' /galaxy-central/config/datatypes_conf.xml \
&& apt-get clean && apt-get purge \
&&  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
ADD topics/tool-generators/docker/welcome.html $TARGDIR/static


ENV GALAXY_CONFIG_BRAND "ToolFactory in Planemo"
EXPOSE 9090

ENTRYPOINT ["/usr/local/bin/planemo" ,"tool_factory", "--galaxy_root" ,"/galaxy-central", "--port", "9090", "--host", "0.0.0.0", "--conda_prefix", "$PDIR/con", "--conda_auto_init", "--conda_auto_install"]

