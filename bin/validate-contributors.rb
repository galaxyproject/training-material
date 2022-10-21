#!/usr/bin/env ruby
require 'yaml'
require 'kwalify'
require './bin/gtn.rb'


# Any error messages
errs = []


CONTRIBUTORS_SCHEMA_UNSAFE = YAML.load_file('bin/schema-contributors.yaml')
CONTRIBUTORS_SCHEMA = automagic_loading(CONTRIBUTORS_SCHEMA_UNSAFE)

# Build validators now that we've filled out the subtopic enum
$contribs_validator = Kwalify::Validator.new(CONTRIBUTORS_SCHEMA)


def validate_document(document, validator)
  errors = validator.validate(document)
  if errors && !errors.empty?
    return errors
  end
  return []
end

errs.push(*validate_document(CONTRIBUTORS, $contribs_validator))

# If we had no errors, validated successfully
if errs.length == 0 then
  puts "\e[38;5;40mCONTRIBUTORS.yaml validated succesfully\e[m"
  exit 0
else
  # Otherwise, print errors and exit non-zero
  puts "\e[48;5;09mCONTRIBUTORS.yaml  has errors\e[m"
  errs.each {|x| puts "  #{x}" }
  exit 1
end

