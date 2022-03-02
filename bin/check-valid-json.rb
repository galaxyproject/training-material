#!/usr/bin/env ruby
require 'json'
fn = ARGV[0]

begin
  JSON.parse(File.open(fn, "r").read)
  puts "#{fn} is valid"
rescue
  puts "#{fn} INVALID"
end
