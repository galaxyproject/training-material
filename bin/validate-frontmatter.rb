#!/usr/bin/env ruby
require 'yaml'
fn = ARGV[0]

# Required keys
tutorial_required_keys = ['layout', 'title', 'time_estimation', 'contributors']
tutorial_optional_keys = ['questions', 'zenodo_link', 'objectives', 'key_points', 'tags', 'edam_ontology', 'requirements', 'follow_up_training']
tutorial_deprecated_keys = ['topic_name', 'tutorial_name', 'type', 'name', 'galaxy_tour', 'hands_on', 'slides', 'workflows']

slides_required_keys = ['layout', 'logo', 'title', 'contributors']
slides_optional_keys = ['time_estimation', 'questions', 'zenodo_link', 'objectives', 'key_points', 'tags', 'edam_ontology', 'requirements', 'follow_up_training', 'class', 'hands_on', 'hands_on_url']
slides_deprecated_keys = ['topic_name', 'tutorial_name', 'type', 'name', 'galaxy_tour', 'slides', 'workflows']

metadata_required_keys = ['name', 'type', 'title', 'summary', 'maintainers']
metadata_optional_keys = ['references', 'requirements', 'docker_image', 'edam_ontology']
metadata_deprecated_keys = ['material']

# Contributors
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')

# Any error messages
errs = []

def skip_disabled(data, fn)
  # If there's an 'enable' key and it is one flavor of 'false', then, exit
  # immediately without testing.
  if data.key?('enable') && (data['enable'] == false || data['enable'].downcase == 'false') then
    puts "#{fn} skipped (disabled)"
    exit 0
  end
end

def check_contributors(data)
  errs = []
  if data.key?('contributors') then
    data['contributors'].each{ |x|
      if not CONTRIBUTORS.key?(x) then
        errs.push("Unknown contributor #{x}, please add to CONTRIBUTORS.yaml")
      end
    }
  end

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

  # Check deprecated keys
  tutorial_deprecated_keys.each{ |x|
    if data.key?(x) then
      errs.push("Deprecated key: #{x}")
    end
  }

  # Check that all keys are valid
  data.keys.each{ |x|
    if not (tutorial_required_keys.include?(x) or tutorial_optional_keys.include?(x) or tutorial_deprecated_keys.include?(x)) then
      errs.push("Unknown key: #{x}")
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
elsif fn.include?('metadata.yaml') then
  data = YAML.load_file(fn)

  # Check for required keys
  metadata_required_keys.each{ |x|
    if not data.key?(x) then
      errs.push("Missing key: #{x}")
    end
  }

  # Check deprecated keys
  metadata_deprecated_keys.each{ |x|
    if data.key?(x) then
      errs.push("Deprecated key: #{x}")
    end
  }

  # Check that all keys are valid
  data.keys.each{ |x|
    if not (metadata_required_keys.include?(x) or metadata_optional_keys.include?(x) or metadata_deprecated_keys.include?(x)) then
      errs.push("Unknown key: #{x}")
    end
  }

elsif fn.include?('slides.html') then
  data = YAML.load_file(fn)
  skip_disabled(data, fn)

  # Check for required keys
  slides_required_keys.each{ |x|
    if not data.key?(x) then
      errs.push("Missing key: #{x}")
    end
  }

  # Check deprecated keys
  slides_deprecated_keys.each{ |x|
    if data.key?(x) then
      errs.push("Deprecated key: #{x}")
    end
  }

  # Check that all keys are valid
  data.keys.each{ |x|
    if not (slides_required_keys.include?(x) or slides_optional_keys.include?(x) or slides_deprecated_keys.include?(x)) then
      errs.push("Unknown key: #{x}")
    end
  }

  # Check that the layout is correct
  if not ['base_slides', 'tutorial_slides'].include?(data['layout']) then
    errs.push("layout should be 'base_slides', not '#{data['layout']}'")
  end

  # Check contributors
  errs = errs.concat(check_contributors(data))
else
  #puts "No validation available for filetype"
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
