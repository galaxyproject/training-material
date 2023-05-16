sample_sentence='Welcome to the Galaxy Training Network! We are a platform for FAIR and Accessible Bioinformatics Training that is open to the world. Most of our content is licensed CC by, 4.0. All of our materials are machine readable through JSON-LD annotations.'
sample=$(mktemp)
echo "$sample_sentence" > "$sample"

ruby bin/ari-synthesize.rb --aws --voice Amy --lang en-GB  --file $sample -o amy.mp3
ruby bin/ari-synthesize.rb --aws --voice Aria --lang en-NZ  --file $sample -o aria.mp3
ruby bin/ari-synthesize.rb --aws --voice Brian --lang en-GB  --file $sample -o brian.mp3
ruby bin/ari-synthesize.rb --aws --voice Emma --lang en-GB  --file $sample -o emma.mp3
ruby bin/ari-synthesize.rb --aws --voice Joanna --lang en-US  --file $sample -o joanna.mp3
ruby bin/ari-synthesize.rb --aws --voice Joey --lang en-US  --file $sample -o joey.mp3
ruby bin/ari-synthesize.rb --aws --voice Kendra --lang en-US  --file $sample -o kendra.mp3
ruby bin/ari-synthesize.rb --aws --voice Matthew --lang en-US  --file $sample -o matthew.mp3
ruby bin/ari-synthesize.rb --aws --voice Nicole --lang en-AU --non-neural --file $sample -o nicole.mp3
ruby bin/ari-synthesize.rb --aws --voice Olivia --lang en-AU  --file $sample -o olivia.mp3
ruby bin/ari-synthesize.rb --aws --voice Raveena --lang en-IN --non-neural --file $sample -o raveena.mp3
ruby bin/ari-synthesize.rb --aws --voice Salli --lang en-US  --file $sample -o salli.mp3
