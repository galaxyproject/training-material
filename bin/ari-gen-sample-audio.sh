sample_sentence='Welcome to the Galaxy Training Network! We are a platform for FAIR and Accessible Bioinformatics Training that is open to the world. Most of our content is licensed CC by, 4.0. All of our materials are machine readable through JSON-LD annotations.'
sample=$(mktemp)
echo "$sample_sentence" > "$sample"

ruby bin/ari-synthesize.rb --aws --voice Amy --lang en-GB  --file $sample -o assets/audio/amy.mp3
ruby bin/ari-synthesize.rb --aws --voice Aria --lang en-NZ  --file $sample -o assets/audio/aria.mp3
ruby bin/ari-synthesize.rb --aws --voice Brian --lang en-GB  --file $sample -o assets/audio/brian.mp3
ruby bin/ari-synthesize.rb --aws --voice Emma --lang en-GB  --file $sample -o assets/audio/emma.mp3
ruby bin/ari-synthesize.rb --aws --voice Joanna --lang en-US  --file $sample -o assets/audio/joanna.mp3
ruby bin/ari-synthesize.rb --aws --voice Joey --lang en-US  --file $sample -o assets/audio/joey.mp3
ruby bin/ari-synthesize.rb --aws --voice Kendra --lang en-US  --file $sample -o assets/audio/kendra.mp3
ruby bin/ari-synthesize.rb --aws --voice Matthew --lang en-US  --file $sample -o assets/audio/matthew.mp3
ruby bin/ari-synthesize.rb --aws --voice Nicole --lang en-AU --non-neural --file $sample -o assets/audio/nicole.mp3
ruby bin/ari-synthesize.rb --aws --voice Olivia --lang en-AU  --file $sample -o assets/audio/olivia.mp3
ruby bin/ari-synthesize.rb --aws --voice Raveena --lang en-IN --non-neural --file $sample -o assets/audio/raveena.mp3
ruby bin/ari-synthesize.rb --aws --voice Salli --lang en-US  --file $sample -o assets/audio/salli.mp3
ruby bin/ari-synthesize.rb --aws --voice Ayanda --lang en-ZA --file $sample -o assets/audio/ayanda.mp3
ruby bin/ari-synthesize.rb --aws --voice Geraint --lang en-GB-WLS --non-neural --file $sample -o assets/audio/geraint.mp3
