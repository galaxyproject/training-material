#!/usr/bin/env ruby
# frozen_string_literal: true

require 'json'

ec = 0
ARGV.each do |fn|
  d = JSON.parse(File.read(fn))

  if !d.is_a?(Hash)
    puts "#{fn}:0:0:e: This notebook is invalid"
    ec = 1
  end

  if d['cells'].empty?
    puts "#{fn}:0:0:e: This notebook is empty"
    ec = 1
  end
rescue StandardError
  puts "#{fn}:0:0:e: This notebook is invalid"
  ec = 1
end

exit ec
