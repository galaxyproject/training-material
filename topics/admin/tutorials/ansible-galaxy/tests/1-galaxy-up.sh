# assumes letsencript staging certificate is preconfigured on VM
galaxy-wait --timeout 30 -g https://$(hostname -f) -v
