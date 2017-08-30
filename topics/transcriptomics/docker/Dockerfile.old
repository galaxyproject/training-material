# Galaxy - RNA-Seq training material

FROM bgruening/galaxy-rna-seq:dev

MAINTAINER Björn A. Grüning, bjoern.gruening@gmail.com

ENV GALAXY_CONFIG_BRAND RNA-Seq training

ADD data_manager_rna_seq.yaml data_manager_rna_seq.yaml

RUN startup_lite && \
    sleep 60 && \
    run-data-managers -u admin@galaxy.org -p admin --config data_manager_rna_seq.yaml -g http://localhost:8080

# Add a Galaxy tour
# ADD <URL to dedicated tour in https://github.com/galaxyproject/galaxy-tours> $GALAXY_ROOT/config/plugins/tours/rna-seq.yaml
