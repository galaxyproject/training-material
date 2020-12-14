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

def synthesize(line, engine)
  digest = Digest::MD5.hexdigest line
  mp3 = File.join(GTN_CACHE, "#{engine}-#{digest}.mp3")
  json = File.join(GTN_CACHE, "#{engine}-#{digest}.json")
  if File.file?(mp3)
    duration = JSON.parse(File.open(json, 'r').read)['end']
    return mp3, json, duration.to_f
  end


  if engine == "aws" then
    # Synthesize
    system('aws', 'polly', 'synthesize-speech', *AWS_PARAMS, '--output-format', 'mp3', '--text', line, mp3) or raise "Failed to speak #{mp3}"
  elsif engine == "mozilla" then
    raw = Tempfile.new('synth-raw')
    system('curl', '--silent', '-G', '--output', raw.path, 'http://localhost:5002/api/tts?text=' + CGI.escape(line))
    system('ffmpeg', '-loglevel', 'error', '-i', raw.path, '-y', mp3)
  end

  # Now collect metadata for JSON
  json_handle = File.open(json, 'w')
  stdout, status = Open3.capture2('ffprobe', '-loglevel', 'warning', '-show_format', '-show_streams', '-print_format', 'json', '-i', mp3)
  data = JSON.parse(stdout)
  duration = data['format']['duration'].to_f

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
