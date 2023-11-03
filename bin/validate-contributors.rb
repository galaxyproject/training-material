#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'kwalify'
require './bin/gtn'

# Any error messages

CONTRIBUTORS_SCHEMA_UNSAFE = YAML.load_file('bin/schema-contributors.yaml')
CONTRIBUTORS_SCHEMA = automagic_loading(CONTRIBUTORS_SCHEMA_UNSAFE)
contribs_validator = Kwalify::Validator.new(CONTRIBUTORS_SCHEMA)

FUNDERS_SCHEMA_UNSAFE = YAML.load_file('bin/schema-funders.yaml')
FUNDERS_SCHEMA = automagic_loading(FUNDERS_SCHEMA_UNSAFE)
funders_validator = Kwalify::Validator.new(FUNDERS_SCHEMA)

ORGANISATIONS_SCHEMA_UNSAFE = YAML.load_file('bin/schema-organisations.yaml')
ORGANISATIONS_SCHEMA = automagic_loading(ORGANISATIONS_SCHEMA_UNSAFE)
organisations_validator = Kwalify::Validator.new(ORGANISATIONS_SCHEMA)

def validate_document(document, validator)
  errors = validator.validate(document)
  return errors if errors && !errors.empty?

  []
end

def show_errors(file, errs)
  # If we had no errors, validated successfully
  if errs.empty?
    puts "\e[38;5;40m#{file} validated succesfully\e[m"
    0
  else
    # Otherwise, print errors and exit non-zero
    puts "\e[48;5;09m#{file}  has errors\e[m"
    errs.each { |x| puts "  #{x}" }
    1
  end
end

ec = 0
# This variable from bin/gtn.rb
errs = validate_document(CONTRIBUTORS, contribs_validator)
ec |= show_errors('CONTRIBUTORS.yaml', errs)
errs = validate_document(FUNDERS, funders_validator)
ec |= show_errors('FUNDERS.yaml', errs)
errs = validate_document(ORGANISATIONS, organisations_validator)
ec | show_errors('ORGANISATIONS.yaml', errs)

# Exit
exit ec
