# Galaxy - Exome sequencing training material

FROM quay.io/bgruening/galaxy:dev

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND Exome sequencing
ENV ENABLE_TTS_INSTALL True

ADD exome_seq.yaml $GALAXY_ROOT/tools.yaml
RUN install-tools $GALAXY_ROOT/tools.yaml
