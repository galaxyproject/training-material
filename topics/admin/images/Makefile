all: ansible-galaxy-0.png ansible-galaxy-1.png  ansible-galaxy-2.png  ansible-galaxy-3.png ansible-galaxy-4.png ansible-galaxy-5.png ansible-galaxy-6.png ansible-galaxy-7.png ansible-galaxy-usegalaxy.png

EXPORT_ARGS := --export-type=png --export-dpi=80 --export-area-page 

clean:
	rm -f ansible-galaxy-0.png ansible-galaxy-1.png  ansible-galaxy-2.png  ansible-galaxy-3.png ansible-galaxy-4.png ansible-galaxy-5.png ansible-galaxy-6.png ansible-galaxy-7.png ansible-galaxy-usegalaxy.png
.PHONY: clean

ansible-galaxy-0.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select ansible-galaxy-1,ansible-galaxy-2,ansible-galaxy-3,cvmfs,connect-to-compute-cluster,pulsar,backup --actions='selection-hide'

ansible-galaxy-1.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select ansible-galaxy-2,ansible-galaxy-3,cvmfs,connect-to-compute-cluster,pulsar,backup --actions='selection-hide'

ansible-galaxy-2.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select ansible-galaxy-3,cvmfs,connect-to-compute-cluster,pulsar,backup --actions='selection-hide'

ansible-galaxy-3.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select cvmfs,connect-to-compute-cluster,pulsar,backup --actions='selection-hide'

ansible-galaxy-4.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select cvmfs,connect-to-compute-cluster,pulsar --actions='selection-hide'

ansible-galaxy-5.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select connect-to-compute-cluster,pulsar --actions='selection-hide'

ansible-galaxy-6.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --select pulsar --actions='selection-hide'

ansible-galaxy-7.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@

ansible-galaxy-usegalaxy.png: schema.svg
	inkscape $(EXPORT_ARGS) schema.svg --export-filename=$@ --actions='select-by-id:minimal;selection-hide;select-by-id:usegalaxy;selection-unhide;'

