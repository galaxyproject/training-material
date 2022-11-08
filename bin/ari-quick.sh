#!/bin/bash
input="$1";
topic=$(echo $input | cut -d/ -f 2)
tutor=$(echo $input | cut -d/ -f 4)
slideslang=$(basename $(echo $input | cut -d/ -f 5) .html)

make _site/training-material/topics/$topic/tutorials/$tutor/$slideslang.pdf ACTIVATE_ENV=pwd

./bin/ari.sh _site/training-material/topics/$topic/tutorials/$tutor/$slideslang.pdf topics/$topic/tutorials/$tutor/$slideslang.html videos/topics/$topic/tutorials/$tutor/$slideslang.mp4
