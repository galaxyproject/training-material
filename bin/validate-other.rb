#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'kwalify'
require './bin/gtn'

# Any error messages
ec = 0
learning_pathway_SCHEMA_UNSAFE = YAML.load_file('bin/schema-learning-pathway.yaml')
learning_pathway_SCHEMA = automagic_loading(learning_pathway_SCHEMA_UNSAFE)

begin
  event_SCHEMA_UNSAFE = YAML.load_file('bin/schema-event.yaml', permitted_classes: [Date])
rescue StandardError
  event_SCHEMA_UNSAFE = YAML.load_file('bin/schema-event.yaml')
end
event_SCHEMA = automagic_loading(event_SCHEMA_UNSAFE)

begin
  event_external_SCHEMA_UNSAFE = YAML.load_file('bin/schema-event-external.yaml', permitted_classes: [Date])
rescue StandardError
  event_external_SCHEMA_UNSAFE = YAML.load_file('bin/schema-event-external.yaml')
end
event_external_SCHEMA = automagic_loading(event_external_SCHEMA_UNSAFE)

Dir.glob('learning-pathways/*.md').reject { |p| p.match(/index.md/) }.each do |file|
  errs = []
  pathway = YAML.load_file(file)

  # Build validators now that we've filled out the subtopic enum
  contribs_validator = Kwalify::Validator.new(learning_pathway_SCHEMA)

  def validate_document(document, validator)
    errors = validator.validate(document)
    return errors if errors && !errors.empty?

    []
  end

  errs.push(*validate_document(pathway, contribs_validator))

  # If we had no errors, validated successfully
  if errs.empty?
    puts "\e[38;5;40m#{file} validated succesfully\e[m"
  else
    # Otherwise, print errors and exit non-zero
    puts "\e[48;5;09m#{file}  has errors\e[m"
    errs.each { |x| puts "  #{x}" }
    ec = 1
  end
end

Dir.glob('events/*.md').reject { |p| p.match(/index.md/) }.each do |file|
  errs = []
  begin
    pathway = YAML.load_file(file, permitted_classes: [Date])
  rescue StandardError
    pathway = YAML.load_file(file)
  end

  # Build validators now that we've filled out the subtopic enum
  contribs_validator = if pathway['layout'] == 'event'
                         Kwalify::Validator.new(event_SCHEMA)
                       else
                         Kwalify::Validator.new(event_external_SCHEMA)
                       end

  def validate_document(document, validator)
    errors = validator.validate(document)
    return errors if errors && !errors.empty?

    []
  end

  errs.push(*validate_document(pathway, contribs_validator))

  # If we had no errors, validated successfully
  if errs.empty?
    puts "\e[38;5;40m#{file} validated succesfully\e[m"
  else
    # Otherwise, print errors and exit non-zero
    puts "\e[48;5;09m#{file}  has errors\e[m"
    errs.each { |x| puts "  #{x}" }
  end
  exit ec
end
