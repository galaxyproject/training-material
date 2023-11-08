#!/usr/bin/env ruby
# frozen_string_literal: true

require 'cgi'
require 'json'
require 'yaml'
require 'optparse'
require 'fileutils'
require 'open3'
require 'tempfile'
require 'digest'
require 'tmpdir'

PUNCTUATION = ['-', '--', '@', '%', '‘', '’', ',', '!', '(', ')', '.', "'", '"', '[', ']', ';', ':'].freeze
ARI_MAP = File.expand_path(File.join(__dir__, 'ari-map.yml'))
WORD_MAP = {}
YAML.load_file(ARI_MAP).each_pair do |k, v|
  WORD_MAP.merge!({ k.downcase => v })
end

GTN_CACHE = File.expand_path(File.join(File.expand_path(__dir__), '..', '.jekyll-cache', 'speech'))
FileUtils.mkdir_p GTN_CACHE

def translate(word)
  return word if /^\s+$/.match(word)

  return word if PUNCTUATION.find_index(word)

  return WORD_MAP[word] if WORD_MAP.key?(word)

  m = /([^A-Za-z0-9]*)([A-Za-z0-9]+)([^A-Za-z0-9]*)(.*)/.match(word)

  if !m
    puts "Error: #{word}"
    return word
  end

  fixed = if m[2]
            WORD_MAP.fetch(m[2].downcase, m[2])
          else
            m[2]
          end

  # puts "#{m} ⇒ #{m[1] + fixed + m[3]}"
  m[1] + fixed + m[3] + m[4]
end

def correct(uncorrected_line)
  # First we try and catch the things we can directly replace (esp usegalaxy.*)
  line = uncorrected_line.strip.split.map do |w|
    translate(w)
  end.join(' ')

  # Now we do more fancy replacements
  line.strip.split(/([ ‘’,'".:;!`()])/).reject(&:empty?).compact.map do |w|
    translate(w)
  end.join
end

def call_engine(engine, line, mp3, voice, lang, neural)
  if engine == 'aws'
    awseng = if neural
               'neural'
             else
               'standard'
             end

    # Synthesize
    args = ['aws', 'polly', 'synthesize-speech', '--engine', awseng, '--language-code', lang, '--voice-id', voice,
            '--output-format', 'mp3', '--text', line, mp3]
    _, stderr, err = Open3.capture3(*args)
    if err.exited? && err.exitstatus.positive?
      puts "ERROR: #{stderr}"
      puts "ERROR: #{err}"
      exit 1
    end
  elsif engine == 'mozilla'
    raw = Tempfile.new('synth-raw')
    _, stderr, err = Open3.capture3('curl', '--silent', '-G', '--output', raw.path,
                                    "http://localhost:5002/api/tts?text=#{CGI.escape(line)}")
    if err.exited? && err.exitstatus.positive?
      puts "ERROR: #{stderr}"
      exit 1
    end

    _, stderr, err = Open3.capture3('ffmpeg', '-loglevel', 'error', '-i', raw.path, '-y', mp3)
    if err.exited? && err.exitstatus.positive?
      puts "ERROR: #{stderr}"
      exit 1
    end
  end
end

def find_duration(mp3)
  stdout, = Open3.capture2('ffprobe', '-loglevel', 'error', '-show_format', '-show_streams', '-print_format', 'json',
                           '-i', mp3)
  data = JSON.parse(stdout)
  data['format']['duration'].to_f
end

def synthesize(uncorrected_line, engine, voice: 'Amy', lang: 'en-GB', neural: true, output: nil)
  line = correct(uncorrected_line)
  digest = Digest::MD5.hexdigest line
  if output.nil?
    mp3 = File.join(GTN_CACHE, "#{engine}-#{digest}-#{voice}.mp3")
    json = File.join(GTN_CACHE, "#{engine}-#{digest}-#{voice}.json")
    if File.file?(mp3)
      duration = JSON.parse(File.read(json))['end']
      return mp3, json, duration.to_f
    end
  else
    mp3 = output
    json = "#{output}.json"
    if File.file?(output)
      return mp3, json, 0.0 # Todo
    end
  end

  # Call our engine
  call_engine(engine, line, mp3, voice, lang, neural)
  duration = find_duration(mp3)

  if line.length < 200 && duration > 27
    # Helena managed to find a specific bad string which, when fed to Mozilla's
    # TTS would generate
    #
    # In: Some important terms you should know.
    # Out Some important terms you should know know know know know know know know know know know know know know ...
    #
    # So we put in a check that the duration hasn't done something crazy, and
    # if it is add something to the end which seems to short-circuit that
    # error.
    #
    # I've reported this upstream but the response was not useful, apparently
    # this is an "expected failure mode".
    #
    # https://github.com/synesthesiam/docker-mozillatts/issues/9
    # https://discourse.mozilla.org/t/sentences-which-trigger-an-endless-loop/72261/8
    warn 'Strange: line was too long'
    call_engine(engine, "#{line}.", mp3)
    duration = find_duration(mp3)
  end

  if line.length < 200 && duration > 27
    # Or maybe they just wrote a super long sentence. Or maybe we need to update the cutoff time.
    warn "ERROR: #{duration} of line is bad: #{line}"
  end

  # Now collect metadata for JSON
  json_handle = File.open(json, 'w')
  json_handle.write(JSON.generate({ time: 0, type: 'sentence', start: 0, end: duration, value: line }))
  json_handle.close

  [mp3, json, duration]
end

def parseOptions
  options = {}
  OptionParser.new do |opts|
    opts.banner = 'Usage: ari-synthesize.rb [options]'

    options[:neural] = true
    options[:voice] = 'Amy'
    options[:lang] = 'en-GB'

    opts.on('--aws', 'Use AWS Polly') do |v|
      options[:aws] = v
    end

    opts.on('--mozilla', 'Use MozillaTTS') do |v|
      options[:mozilla] = v
    end

    opts.on('--non-neural', '[AWS] Non-neural voice') do |_v|
      options[:neural] = false
    end

    opts.on('--voice=VOICE', '[AWS] Voice ID') do |n|
      options[:voice] = n
    end

    opts.on('--lang=LANG', '[AWS] Language code') do |n|
      options[:lang] = n
    end

    opts.on('-fFILE', '--file=FILE', 'File containing line of text to speak') do |n|
      options[:file] = n
    end

    opts.on('-oFILE', '--output=FILE', 'Location to save the file in (defaults to auto-generated location)') do |n|
      options[:output] = n
    end

    opts.on('-v', '--[no-]verbose', 'Run verbosely') do |v|
      options[:verbose] = v
    end
  end.parse!

  if !(options[:aws] || options[:mozilla])
    puts 'ERROR: You must use aws or mozilla'
    exit 1
  end

  if !(options[:file])
    puts 'ERROR: You must provide a file with a single sentence to speak'
    exit 1
  end

  sentence = File.read(options[:file]).chomp
  if options[:aws]
    engine = 'aws'
  elsif options[:mozilla]
    engine = 'mozilla'
  end

  [sentence, engine, options]
end

if __FILE__ == $PROGRAM_NAME
  sentence, engine, options = parseOptions
  mp3, = synthesize(sentence, engine, voice: options[:voice], lang: options[:lang], neural: options[:neural],
                                      output: options[:output])
  puts mp3
end
