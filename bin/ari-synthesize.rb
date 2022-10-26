#!/usr/bin/env ruby
require 'cgi'
require 'json'
require 'yaml'
require 'optparse'
require 'fileutils'
require 'open3'
require 'tempfile'
require 'digest'
require 'tmpdir'

PUNCTUATION = ['-', '--', '@', '%', '‘', '’', ',', '!', '(', ')', '.', "'", '"', '[', ']', ';', ':']
ARI_MAP = File.expand_path(File.join(__dir__, 'ari-map.yml'))
WORD_MAP = {}
YAML.load_file(ARI_MAP).each_pair do |k,v|
 WORD_MAP.merge!({k.downcase => v})
end

GTN_CACHE = File.expand_path(File.join(File.expand_path(__dir__), '..', '.jekyll-cache', 'speech'))
FileUtils.mkdir_p GTN_CACHE

def translate(word)
  if /^\s+$/.match(word)
    return word
  end

  if PUNCTUATION.find_index(word) then
    return word
  end

  if WORD_MAP.key?(word) then
    return WORD_MAP[word]
  end

  m = /([^A-Za-z0-9]*)([A-Za-z0-9]+)([^A-Za-z0-9]*)(.*)/.match(word)

  if ! m then
    puts "Error: #{word}"
    return word
  end

  if m[2] then
    fixed = WORD_MAP.fetch(m[2].downcase, m[2])
  else
    fixed = m[2]
  end

  #puts "#{m} ⇒ #{m[1] + fixed + m[3]}"
  return m[1] + fixed + m[3] + m[4]
end

def correct(uncorrected_line)
  # First we try and catch the things we can directly replace (esp usegalaxy.*)
  line = uncorrected_line.strip.split(' ').map{ |w|
    translate(w)
  }.join(' ')

  # Now we do more fancy replacements
  line = line.strip.split(/([ ‘’,'".:;!`()])/).reject(&:empty?).compact.map{ |w|
    translate(w)
  }.join('')
  line
end

def call_engine(engine, line, mp3, voice, lang, neural)
  if engine == "aws" then
    if neural then
      awseng = 'neural'
    else
      awseng = 'standard'
    end

    # Synthesize
    args = ['aws', 'polly', 'synthesize-speech', "--engine", awseng, "--language-code", lang, "--voice-id", voice, '--output-format', 'mp3', '--text', line, mp3]
    _, stderr, err = Open3.capture3(*args)
    if err.exited? and err.exitstatus > 0
      puts "ERROR: #{stderr}"
      puts "ERROR: #{err}"
      exit 1
    end
  elsif engine == "mozilla" then
    raw = Tempfile.new('synth-raw')
    _, stderr, err = Open3.capture3('curl', '--silent', '-G', '--output', raw.path, 'http://localhost:5002/api/tts?text=' + CGI.escape(line))
    if err.exited? and err.exitstatus > 0
      puts "ERROR: #{stderr}"
      exit 1
    end

    _, stderr, err = Open3.capture3('ffmpeg', '-loglevel', 'error', '-i', raw.path, '-y', mp3)
    if err.exited? and err.exitstatus > 0
      puts "ERROR: #{stderr}"
      exit 1
    end
  end
end

def find_duration(mp3)
  stdout, _ = Open3.capture2('ffprobe', '-loglevel', 'error', '-show_format', '-show_streams', '-print_format', 'json', '-i', mp3)
  data = JSON.parse(stdout)
  duration = data['format']['duration'].to_f
  return duration
end

def synthesize(uncorrected_line, engine, voice: "Amy", lang: "en-GB", neural: true)
  line = correct(uncorrected_line)
  digest = Digest::MD5.hexdigest line
  mp3 = File.join(GTN_CACHE, "#{engine}-#{digest}-#{voice}.mp3")
  json = File.join(GTN_CACHE, "#{engine}-#{digest}-#{voice}.json")
  if File.file?(mp3)
    duration = JSON.parse(File.open(json, 'r').read)['end']
    return mp3, json, duration.to_f
  end

  # Call our engine
  call_engine(engine, line, mp3, voice, lang, neural)
  duration = find_duration(mp3)

  if line.length < 200 && duration > 27
    # Helena managed to find a specific bad string which, when fed to Mozilla's
    # TTS would generate
    #
    # In: Some important terms you should know.
    # Out Some important terms you should know know know know know know know know know know know know know know know know ...
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
    STDERR.puts "Strange: line was too long"
    call_engine(engine, line + '.', mp3)
    duration = find_duration(mp3)
  end

  if line.length < 200 && duration > 27
    # Or maybe they just wrote a super long sentence. Or maybe we need to update the cutoff time.
    STDERR.puts "ERROR: #{duration} of line is bad: #{line}"
  end

  # Now collect metadata for JSON
  json_handle = File.open(json, 'w')
  json_handle.write(JSON.generate({"time": 0, "type": "sentence", "start": 0, "end": duration, "value": line}))
  json_handle.close

  return mp3, json, duration
end

def parseOptions
  options = {}
  OptionParser.new do |opts|
    opts.banner = "Usage: ari-synthesize.rb [options]"

    options[:neural] = true
    options[:voice] = 'Amy'
    options[:lang] = 'en-GB'

    opts.on("--aws", "Use AWS Polly") do |v|
      options[:aws] = v
    end

    opts.on("--mozilla", "Use MozillaTTS") do |v|
      options[:mozilla] = v
    end

    opts.on("--non-neural", "[AWS] Non-neural voice") do |v|
      options[:neural] = false
    end

    opts.on("--voice=VOICE", "[AWS] Voice ID") do |n|
      options[:voice] = n
    end

    opts.on("--lang=LANG", "[AWS] Language code") do |n|
      options[:lang] = n
    end

    opts.on("-fFILE", "--file=FILE", "File containing line of text to speak") do |n|
      options[:file] = n
    end

    opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
      options[:verbose] = v
    end
  end.parse!

  if not (options[:aws] or options[:mozilla])
    puts "ERROR: You must use aws or mozilla"
    exit 1
  end

  if not options[:file]
    puts "ERROR: You must provide a file with a single sentence to speak"
    exit 1
  end

  sentence = File.open(options[:file], 'r').read.chomp
  if options[:aws]
    engine = 'aws'
  elsif options[:mozilla]
    engine = 'mozilla'
  end

  return sentence, engine, options
end

if __FILE__ == $0
  sentence, engine, options = parseOptions()
  mp3, _ = synthesize(sentence, engine, voice: options[:voice], lang: options[:lang], neural: options[:neural])
  puts mp3
end
