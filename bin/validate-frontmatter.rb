#!/usr/bin/env ruby
require 'yaml'
require 'pathname'
require 'kwalify'
fn = ARGV[0]

metadata_schema = YAML.load_file('bin/schema-topic.yaml')
tutorial_schema = YAML.load_file('bin/schema-tutorial.yaml')
slides_schema = YAML.load_file('bin/schema-slides.yaml')

# Contributors
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')

# Any error messages
errs = []

data = YAML.load_file(fn)

# If it's disabled, exit early
if data.key?('enable') && (data['enable'] == false || data['enable'].downcase == 'false') then
  puts "#{fn} skipped (disabled)"
  exit 0
end

if not fn.include?('metadata.yaml') then
  # Load topic metadata
  topic = fn.split('/')[1]
  topic_metadata = YAML.load_file("topics/#{topic}/metadata.yaml")

  # Load subtopic titles
  if data.key?('subtopic') then
    subtopic_ids = []
    topic_metadata['subtopics'].each{ |x|
      subtopic_ids.push(x['id'])
    }

    tutorial_schema['mapping']['subtopic']['enum'] = subtopic_ids
    slides_schema['mapping']['subtopic']['enum'] = subtopic_ids
  end
end

# Build validators now that we've filled out the subtopic enum
metadata_validator = Kwalify::Validator.new(metadata_schema)
tutorial_validator = Kwalify::Validator.new(tutorial_schema)
slides_validator = Kwalify::Validator.new(slides_schema)

def check_contributors(input)
  errs = []
  if input.key?('contributors') then
    key = 'contributors'
  elsif input.key?('maintainers') then
    key = 'maintainers'
  end

  input[key].each{ |x|
    if not CONTRIBUTORS.key?(x) then
      errs.push("Unknown #{key} #{x}, please add to CONTRIBUTORS.yaml")
    end
  }

  return errs
end

def validate_non_empty_key_value(map, key)
    if map.key?(key) then
      if map[key].length == 0 then
        return ["Empty #{key} for requirement"]
      end
    else
      return ["Missing #{key} for requirement"]
    end
    return []
end

def validate_requirements(requirements)
  errs = []
  # Exit early if no requirements
  if requirements.nil? or requirements.length == 0
    return []
  end

  # Otherwise check each
  for requirement in requirements
    # For external links, they need a link that is non-empty
    if requirement['type'] == 'external'
      errs.push(*validate_non_empty_key_value(requirement, 'title'))
      errs.push(*validate_non_empty_key_value(requirement, 'link'))

      requirement.keys.each{ |x|
        if not ['title', 'link', 'type'].include?(x) then
          errs.push("Unknown key #{x}")
        end
      }
    elsif requirement['type'] == 'internal'
      errs.push(*validate_non_empty_key_value(requirement, 'topic_name'))
      errs.push(*validate_non_empty_key_value(requirement, 'tutorials'))

      requirement.keys.each{ |x|
        if not ['topic_name', 'tutorials', 'type'].include?(x) then
          errs.push("Unknown key #{x}")
        end
      }
      # For the internal requirements, test that they point at something real.
      if requirement.key?('tutorials') then
        requirement['tutorials'].each{ |tutorial|
          # For each listed tutorial check that a directory with that name exists
          pn = Pathname.new("topics/#{requirement['topic_name']}/tutorials/#{tutorial}")

          if not pn.directory?
            errs.push("Internal requirement to topics/#{requirement['topic_name']}/tutorials/#{tutorial} does not exist")
          end
        }
      end
      #
    elsif requirement['type'] == 'none'
      errs.push(*validate_non_empty_key_value(requirement, 'title'))

      requirement.keys.each{ |x|
        if not ['title', 'type'].include?(x) then
          errs.push("Unknown key #{x}")
        end
      }
    else
      errs.push("Unknown requirement type #{requirement['type']}")
    end
  end

  return errs
end

# Handle tutorials
if fn.include?('tutorial.md') then
  # Validate document
  errors = tutorial_validator.validate(data)
  if errors && !errors.empty?
    for e in errors
      errs.push("[#{e.path}] #{e.message}")
    end
  end

  # Check requirements
  if data.key?('requirements') then
    errs.push(*validate_requirements(data['requirements']))
  end

  # Check follow ups
  if data.key?('follow_up_training') then
    errs.push(*validate_requirements(data['follow_up_training']))
  end

  # Check contributors
  errs = errs.concat(check_contributors(data))


# Validate Metadata
elsif fn.include?('metadata.yaml') then
  # Validate document
  errors = metadata_validator.validate(data)
  if errors && !errors.empty?
    for e in errors
      errs.push("[#{e.path}] #{e.message}")
    end
  end

  # Check contributors
  errs = errs.concat(check_contributors(data))

elsif fn.include?('slides.html') then
  # Validate document
  errors = slides_validator.validate(data)
  if errors && !errors.empty?
    for e in errors
      errs.push("[#{e.path}] #{e.message}")
    end
  end

  # Check requirements
  if data.key?('requirements') then
    errs.push(*validate_requirements(data['requirements']))
  end

  # Check contributors
  errs = errs.concat(check_contributors(data))
else
  #puts "No validation available for filetype"
  exit 0
end


# If we had no errors, validated successfully
if errs.length == 0 then
  puts "\e[38;5;40m#{fn} validated succesfully\e[m"
  exit 0
else
  # Otherwise, print errors and exit non-zero
  puts "\e[48;5;09m#{fn} has errors\e[m"
  errs.each {|x| puts "  #{x}" }
  exit 1
end
