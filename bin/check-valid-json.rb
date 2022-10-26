#!/usr/bin/env ruby
require 'json'

ARGV.each{|fn|
  begin
    JSON.parse(File.open(fn, "r").read)
  rescue
    puts "#{fn}:0:0:e: This notebook is invalid"
  end
}
