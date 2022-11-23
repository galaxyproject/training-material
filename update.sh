#update the copy served on ross's home server 
sudo cp -r /home/ross/rossgit/training-material/_site/* /var/www/html/
sudo sed -iv 's+href="/training-material/videos/watch.html?v=/galaxy-project/slides/introduction"+href="/training-material/videos/introduction.mp4"+g' \
 /var/www/html/training-material/topics/galaxy-project/index.html
 
sudo sed -iv 's+href="/training-material/videos/watch.html?v=/galaxy-project/slides/introduction"+href="/training-material/videos/introduction.mp4"+g' \
 /var/www/html/training-material/topics/galaxy-project/slides/introduction.html
 
sudo chown -R www-data /var/www/html/training-material
