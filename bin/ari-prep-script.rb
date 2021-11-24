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

END_OF_SENTENCE_DURATION = 0.2
END_OF_SLIDE_DURATION = 1.04

editly = {
  'width' => 1920,
  'height' => 1080,
  'fps' => 30,
  'fast' => ENV.fetch('EDITLY_FAST', "false") == "true",
  'outPath' => File.join(dir, 'tmp.mp4'),
  'defaults' => {
    'transition' => {
      'duration': 0,
    },
  },
  'keepSourceAudio': true,
  'clips' => []
}

# Intro slides.
# Fast: editly --json editly.json5  126,23s user 5,62s system 126% cpu 1:44,08 total
# Slow: editly --json editly.json5  902,71s user 69,27s system 326% cpu 4:57,54 total

def timefmt(t, fmt)
  seconds = t % (24 * 3600)
  hours = seconds.to_i / 3600
  seconds = seconds % 3600
  minutes = seconds.to_i / 60
  seconds = seconds % 60
  (seconds, ms) = seconds.divmod(1)
  seconds = seconds
  ms = 1000 * ms

  if fmt == "vtt"
    "%02d:%02d:%02d.%03d" % [hours, minutes, seconds, ms]
  else
    "%02d:%02d:%02d,%03d" % [hours, minutes, seconds, ms]
  end
end

def split_sentence(sentence, timing)
  split_sentence = sentence.split(' ')

  chunks = split_sentence.each_slice(20).to_a.length
  split_sentence.each_slice(20).each_with_index.map{|chunk, idx|
    t0 = timing * (idx / chunks)
    tf = timing * ((1 + idx) / chunks)
   [chunk, t0, tf]
  }
end

# For each line in the script, iterate
# The spoken line is the one that's hashed. (2nd col)
editly['clips'] = script['blocks'].map.with_index{ |phrases, idx|
  # We're inside a single slide.
  parts = phrases.map{|subtitle|
    digest = Digest::MD5.hexdigest subtitle

    # Write out our script bits.
    handle = File.open(File.join(dir, digest + '.txt'), 'w')
    handle.write(subtitle)

    handle = File.open(File.join(dir, digest + '-subtitle.txt'), 'w')
    handle.write(subtitle)

    # Synthesize and copy to the temp dir
    voice = script['voice']

    mp3, json, duration = synthesize(subtitle, engine, 'voice': voice['id'], 'lang': voice['lang'], neural: voice['neural'])
    puts "\tSynthesizing: #{subtitle}"
    FileUtils.cp(mp3, File.join(dir, digest + '.mp3'))
    FileUtils.cp(json, File.join(dir, digest + '.json'))

    [
      {
        :caption => subtitle,
        'duration' => duration,
        'layers' => [{
          'type' => 'audio',
          'path' => File.join(dir, "#{digest}.mp3"),
        }, {
          'type' => 'image',
          'path' => File.join(dir, "slides.%03d.png" % idx),
          'resizeMode' => 'stretch',
          'zoomDirection' => nil,
        }]
      }, {
        'duration' => END_OF_SENTENCE_DURATION,
        'layers' => [{
          'type' => 'image',
          'path' => File.join(dir, "slides.%03d.png" % idx),
          'resizeMode' => 'stretch',
          'zoomDirection' => nil,
        }]
      }
    ]
  }
  parts.flatten!
  # Strip out the last pause
  parts = parts[0..-2]

  # Here we add 1 second of silence at the end of each slide.
  parts.push({
    'transition' => {
      'name' => 'fadegrayscale',
      'duration' => 2,
    },
    'duration' => END_OF_SLIDE_DURATION,
    'layers' => [{
      'type' => 'image',
      'path' => File.join(dir, "slides.%03d.png" % idx),
      'resizeMode' => 'stretch',
      'zoomDirection' => nil,
    }]
  })
  parts.flatten
}.flatten


subtitle_timings = []
offset = 0
editly['clips'].each{|layer|
  if layer.has_key?(:caption)
    subtitle_timings += split_sentence(layer[:caption], layer['duration']).map{|sen_part, time_prev, time_next|
      [sen_part.join(' '), offset + time_prev, offset + time_next]
    }
  end
  offset += layer['duration']
}

# Remove our :caption key.
editly['clips'].map!{|layer|
  layer.delete(:caption)
  layer
}

video_script = File.open(File.join(dir, 'editly.json5'), 'w')
video_script.write(JSON.generate(editly))

vtt = File.open(File.join(dir, 'out.vtt'), 'w')
srt = File.open(File.join(dir, 'out.srt'), 'w')

vtt.write("WEBVTT\n\n\n")
subtitle_timings.each_with_index{|subtitle, index|
  sub, time_prev, time_next = subtitle

  vtt.write("#{index}\n")
  srt.write("#{index}\n")
  vtt.write("#{timefmt(time_prev, 'vtt')} --> #{timefmt(time_next, 'vtt')}\n")
  srt.write("#{timefmt(time_prev, 'srt')} --> #{timefmt(time_next, 'srt')}\n")
  vtt.write("#{sub}\n\n")
  srt.write("#{sub}\n\n")
}
