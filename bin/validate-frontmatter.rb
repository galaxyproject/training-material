#!/usr/bin/env ruby
require 'yaml'
fn = ARGV[0]

# Required keys
tutorial_required_keys = ['layout', 'title', 'time_estimation', 'contributors']
slides_required_keys = ['layout', 'logo', 'title', 'contributors']

# Contributors
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')

# Any error messages
errs = []

def skip_disabled(data, fn)
  # If there's an 'enable' key and it is one flavor of 'false', then, exit
  # immediately without testing.
  if data.key?('enable') && (data['enable'].downcase == 'false' || data['enable'] == false) then
    puts "#{fn} skipped (disabled)"
    exit 0
  end
end

def check_contributors(data)
  errs = []
  data['contributors'].each{ |x|
    if not CONTRIBUTORS.key?(x) then
      errs.push("Unknown contributor #{x}, please add to CONTRIBUTORS.yaml")
    end
  }

  return errs
end

# Handle tutorials
if fn.include?('tutorial.md') then
  data = YAML.load_file(fn)
  skip_disabled(data, fn)

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

  # Check contributors
  errs = errs.concat(check_contributors(data))
elsif fn.include?('slides.html') then
  data = YAML.load_file(fn)
  skip_disabled(data, fn)

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

  # Check contributors
  errs = errs.concat(check_contributors(data))
else
  puts "No validation available for filetype"
  exit 0
end


# If we had no errors, validated successfully
if errs.length == 0 then
  puts "#{fn} validated succesfully"
  exit 0
else
  # Otherwise, print errors and exit non-zero
  puts "#{fn} has errors"
  errs.each {|x| puts "  #{x}" }
  exit 1
end
