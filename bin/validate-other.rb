#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'kwalify'
require './bin/gtn'

# Any error messages

learning_pathway_SCHEMA_UNSAFE = YAML.load_file('bin/schema-learning-pathway.yaml')
learning_pathway_SCHEMA = automagic_loading(learning_pathway_SCHEMA_UNSAFE)

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
  end
end
