#!/usr/bin/env ruby
# Given a script (json file)
# explode it into a directory of individual lines in .txt format (.txt, -sub.txt)
# For each spoken line convert into mp3
# Then re-assemble a script for ffmpeg
require File.expand_path(File.dirname(__FILE__)) + '/ari-synthesize.rb'
require 'fileutils'
require 'digest'
require 'json'

script = JSON.parse(File.open(ARGV[0], 'r').read)
dir = ARGV[1]
engine = ARGV[2]

# For each line in the script, iterate
# The spoken line is the one that's hashed. (2nd col)
final_script = script.map.with_index{ |slide, idx|
  paired = slide[0].zip(slide[1])
  paired = paired.map{ |subtitle, spoken|
    digest = Digest::MD5.hexdigest spoken

    # Write out our script bits.
    handle = File.open(File.join(dir, digest + '.txt'), 'w')
    handle.write(spoken)

    handle = File.open(File.join(dir, digest + '-subtitle.txt'), 'w')
    handle.write(subtitle)

    # Synthesize and copy to the temp dir
    mp3, json, duration = synthesize(spoken, engine)
    puts "\tSynthesizing: #{spoken}"
    FileUtils.cp(mp3, File.join(dir, digest + '.mp3'))
    FileUtils.cp(json, File.join(dir, digest + '.json'))

    [
      "file 'slides.%03d.png'\nduration #{duration}\n" % idx,
      "file '#{digest}.mp3'\nduration #{duration}\n",
    ]
  }
  # Silence at the end of the slide
  paired.push([
      "file 'slides.%03d.png'\nduration 1.04\n" % idx,
      "file 'silence.mp3'\nduration 1.04\n",
  ])

  paired
}

video_script = File.open(File.join(dir, 'images.txt'), 'w')
audio_script = File.open(File.join(dir, 'sounds.txt'), 'w')

video_script.write(
  final_script.map{|slide|
    slide.map{|video, audio|
      video
    }.join('')
  }.join('')
)

audio_script.write(
  final_script.map{|slide|
    slide.map{|video, audio|
      audio
    }.join('')
  }.join('')
)
