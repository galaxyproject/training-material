#!/usr/bin/env ruby
require 'yaml'
fn = ARGV[0]

# Required keys
tutorial_required_keys = ['layout', 'title', 'time_estimation', 'contributors']
slides_required_keys = ['layout', 'logo', 'title', 'contributors']

# Any error messages
errs = []

# Handle tutorials
if fn.include?('tutorial.md') then
  data = YAML.load_file(fn)

  # Check for required keys
  tutorial_required_keys.each{ |x|
    if not data.key?(x) then
      errs.push("Missing key: #{x}")
    end
  }

  # Check that the layout is correct
  if data['layout'] != "tutorial_hands_on" then
    errs.push("layout should be 'tutorial_hands_on', not '#{data['layout']}'")
  end

  # Check time formatting
  if data.key?('time_estimation') then
    match = /^(?:([0-9]*)[Hh])*(?:([0-9]*)[Mm])*(?:([0-9.]*)[Ss])*$/.match(data['time_estimation'])
    if match.nil? then
      errs.push("Time specification could not be parsed (Should be of form ##h##m##s, is '#{data['time_estimation']}')")
    end
  end

elsif fn.include?('slides.html') then
  data = YAML.load_file(fn)

  # Check for required keys
  slides_required_keys.each{ |x|
    if not data.key?(x) then
      errs.push("Missing key: #{x}")
    end
  }

  # Check that the layout is correct
  if not ['base_slides', 'tutorial_slides'].include?(data['layout']) then
    errs.push("layout should be 'base_slides', not '#{data['layout']}'")
  end
else
  puts "No validation available for filetype"
  exit 0
end

if errs.length == 0 then
  puts "#{fn} validated succesfully"
  exit 0
else
  puts "#{fn} has errors"
  errs.each {|x| puts "  #{x}" }
  exit 1
end
