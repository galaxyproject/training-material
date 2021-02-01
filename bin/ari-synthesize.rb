#!/usr/bin/env ruby
require 'cgi'
require 'json'
require 'optparse'
require 'fileutils'
require 'open3'
require 'tempfile'
require 'digest'
require 'tmpdir'

GTN_CACHE = File.expand_path(File.join(File.expand_path(__dir__), '..', '.jekyll-cache', 'speech'))
AWS_PARAMS = ["--engine", "neural", "--language-code", "en-GB", "--voice-id", "Amy"]
FileUtils.mkdir_p GTN_CACHE

def call_engine(engine, line, mp3)
  if engine == "aws" then
    # Synthesize
    _, stderr, err = Open3.capture3('aws', 'polly', 'synthesize-speech', *AWS_PARAMS, '--output-format', 'mp3', '--text', line, mp3)
    if err.exited? and err.exitstatus > 0
      puts "ERROR: #{stderr}"
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

def synthesize(line, engine)
  digest = Digest::MD5.hexdigest line
  mp3 = File.join(GTN_CACHE, "#{engine}-#{digest}.mp3")
  json = File.join(GTN_CACHE, "#{engine}-#{digest}.json")
  if File.file?(mp3)
    duration = JSON.parse(File.open(json, 'r').read)['end']
    return mp3, json, duration.to_f
  end

  # Call our engine
  call_engine(engine, line, mp3)
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

    opts.on("--aws", "Use AWS Polly") do |v|
      options[:aws] = v
    end

    opts.on("--mozilla", "Use MozillaTTS") do |v|
      options[:mozilla] = v
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

  return sentence, engine
end

if __FILE__ == $0
  sentence, engine = parseOptions()
  mp3, _ = synthesize(sentence, engine)
  puts mp3
end
