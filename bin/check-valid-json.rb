#!/usr/bin/env ruby
# frozen_string_literal: true

require 'json'

ARGV.each do |fn|
  JSON.parse(File.read(fn))
rescue StandardError
  puts "#{fn}:0:0:e: This notebook is invalid"
end
