#!/usr/bin/env ruby
require 'json'

ARGV.each do |fn|
  JSON.parse(File.read(fn))
rescue StandardError
  puts "#{fn}:0:0:e: This notebook is invalid"
end
