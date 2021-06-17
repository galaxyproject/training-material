#!/bin/bash
input="$1";
topic=$(echo $input | cut -d/ -f 2)
tutor=$(echo $input | cut -d/ -f 4)

make _site/training-material/topics/$topic/tutorials/$tutor/slides.pdf ACTIVATE_ENV=pwd

./bin/ari.sh _site/training-material/topics/$topic/tutorials/$tutor/slides.pdf topics/$topic/tutorials/$tutor/slides.html videos/topics/$topic/tutorials/$tutor/slides.mp4
