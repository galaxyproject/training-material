# assumes letsencript staging certificate is preconfigured on VM
galaxy-wait --timeout 30 -g https://"${GALAXY_HOSTNAME}" -v
