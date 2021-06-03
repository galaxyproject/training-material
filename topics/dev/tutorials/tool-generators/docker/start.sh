# please adjust as needed and save as start.sh
# the volume can be populated with newly generated tools from unpacked toolshed archives and remains persistent
# history WILL NOT persist so save it if you want to import it in the future to save typing all your work again.
# first run will take ages to fill the local galaxy-central volume for persistence
sudo docker run -d --network host --privileged  \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /home/ross/export:/export \
    quay.io/fubar2/toolfactory_tutorial:latest

