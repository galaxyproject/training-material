---
title: "`docker-compose up` fails with error `/usr/bin/start.sh: line 133: /galaxy/.venv/bin/uwsgi: No such file or directory`"
box_type: question
layout: faq
contributors: [fubar2]
---


- This is why it's useful to watch the boot process without detaching
- This can happen if a container has become corrupt on disk after being interrupted
  - cured by a complete cleanup.
- Make sure no docker galaxy-server related processes are running - use docker ps to check and stop them manually
- delete the `..compose/export directory` with `sudo rm -rf export/*` to clean out any corrupted files
- run `docker system prune` to clear out any old corrupted containers, images or networks. Then run `docker volume prune` in the same way to remove the shared volumes.
- run `docker-compose pull` again to ensure the images are correct
- run `docker-compose up` to completely rebuild the appliance from scratch. **Please be patient**.


