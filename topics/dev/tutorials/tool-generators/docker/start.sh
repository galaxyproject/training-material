# please adjust as needed and save as start.sh
# the volume can be populated with newly generated tools from unpacked toolshed archives and remains persistent
# history WILL NOT persist so save it if you want to import it in the future to save typing all your work again.
docker run -d -p 8080:80 -p 9090:9090  \
    -v /home/ross/rossgit/training-material/mytools:/planemo/mytools \
    tf-tutorial

