#!/usr/bin/env ruby
require 'yaml'
require 'pathname'
require 'kwalify'
fn = ARGV[0]

metadata_schema = YAML.load_file('bin/schema-topic.yaml')
tutorial_schema = YAML.load_file('bin/schema-tutorial.yaml')
slides_schema = YAML.load_file('bin/schema-slides.yaml')
requirement_external_schema = YAML.load_file('bin/schema-requirement-external.yaml')
requirement_internal_schema = YAML.load_file('bin/schema-requirement-internal.yaml')

# Any error messages
errs = []

data = YAML.load_file(fn)

# Contributors
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')
# Update the existing schemas to have enums with values. Then we get validation *for free*!
tutorial_schema['mapping']['contributors']['sequence'][0]['enum'] = CONTRIBUTORS.keys
slides_schema['mapping']['contributors']['sequence'][0]['enum'] = CONTRIBUTORS.keys
metadata_schema['mapping']['maintainers']['sequence'][0]['enum'] = CONTRIBUTORS.keys

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
$metadata_validator = Kwalify::Validator.new(metadata_schema)
$tutorial_validator = Kwalify::Validator.new(tutorial_schema)
$slides_validator = Kwalify::Validator.new(slides_schema)
$requirement_external_validator = Kwalify::Validator.new(requirement_external_schema)
$requirement_internal_validator = Kwalify::Validator.new(requirement_internal_schema)


def validate_document(document, validator)
  errors = validator.validate(document)
  if errors && !errors.empty?
    return errors
  end
  return []
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
      errs.push(*validate_document(requirement, $requirement_external_validator))
    elsif requirement['type'] == 'internal'
      errs.push(*validate_document(requirement, $requirement_internal_validator))

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

# Generic error handling:
## Check requirements
if data.key?('requirements') then
  errs.push(*validate_requirements(data['requirements']))
end

## Check follow ups
if data.key?('follow_up_training') then
  errs.push(*validate_requirements(data['follow_up_training']))
end

# Custom error handling:
if fn.include?('tutorial.md') then
  errs.push(*validate_document(data, $tutorial_validator))
elsif fn.include?('metadata.yaml') then
  errs.push(*validate_document(data, $metadata_validator))
elsif fn.include?('slides.html') then
  errs.push(*validate_document(data, $slides_validator))
else
  errs.push("No validation available for this type of file")
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
