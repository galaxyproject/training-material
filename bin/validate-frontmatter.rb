#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'find'
require 'pathname'
require 'kwalify'
require './bin/gtn'

# Validates the frontmatter of all files
module SchemaValidator
  # Schemas
  @TOPIC_SCHEMA_UNSAFE = YAML.load_file('bin/schema-topic.yaml')
  @TUTORIAL_SCHEMA_UNSAFE = YAML.load_file('bin/schema-tutorial.yaml')
  @SLIDES_SCHEMA_UNSAFE = YAML.load_file('bin/schema-slides.yaml')
  @FAQ_SCHEMA_UNSAFE = YAML.load_file('bin/schema-faq.yaml')
  @QUIZ_SCHEMA_UNSAFE = YAML.load_file('bin/schema-quiz.yaml')
  @requirement_external_schema = YAML.load_file('bin/schema-requirement-external.yaml')
  @requirement_internal_schema = YAML.load_file('bin/schema-requirement-internal.yaml')

  # Update the existing schemas to have enums with values. Then we get validation *for free*!
  @TUTORIAL_SCHEMA = automagic_loading(@TUTORIAL_SCHEMA_UNSAFE)
  @SLIDES_SCHEMA = automagic_loading(@SLIDES_SCHEMA_UNSAFE)
  @TOPIC_SCHEMA = automagic_loading(@TOPIC_SCHEMA_UNSAFE)
  @FAQ_SCHEMA = automagic_loading(@FAQ_SCHEMA_UNSAFE)
  @QUIZ_SCHEMA = automagic_loading(@QUIZ_SCHEMA_UNSAFE)

  @TUTORIAL_SCHEMA['mapping']['contributions']['required'] = false
  @SLIDES_SCHEMA['mapping']['contributions']['required'] = false

  # Build validators now that we've filled out the subtopic enum
  @topic_validator = Kwalify::Validator.new(@TOPIC_SCHEMA)
  @tutorial_validator = Kwalify::Validator.new(@TUTORIAL_SCHEMA)
  @slides_validator = Kwalify::Validator.new(@SLIDES_SCHEMA)
  @faq_validator = Kwalify::Validator.new(@FAQ_SCHEMA)
  @quiz_validator = Kwalify::Validator.new(@QUIZ_SCHEMA)
  @requirement_external_validator = Kwalify::Validator.new(@requirement_external_schema)
  @requirement_internal_validator = Kwalify::Validator.new(@requirement_internal_schema)

  def self.validate_document(document, validator)
    errors = validator.validate(document)
    return errors if errors && !errors.empty?

    []
  end

  def self.validate_non_empty_key_value(map, key)
    return ["Missing #{key} for requirement"] unless map.key?(key)
    return ["Empty #{key} for requirement"] if map[key].empty?

    []
  end

  def self.tutorial?(fn)
    fn.include?('tutorial.md') || fn =~ /tutorial_[A-Z]{2,}.md/
  end

  def self.slide?(fn)
    fn.include?('slides.html') || fn =~ /slides_[A-Z]{2,}.html/
  end

  def self.validate_requirements(requirements)
    errs = []
    # Exit early if no requirements
    return [] if requirements.nil? || requirements.empty?

    # Otherwise check each
    requirements.each do |requirement|
      # For external links, they need a link that is non-empty
      case requirement['type']
      when 'external'
        errs.push(*validate_document(requirement, @requirement_external_validator))
      when 'internal'
        errs.push(*validate_document(requirement, @requirement_internal_validator))

        # For the internal requirements, test that they point at something real.
        if requirement.key?('tutorials')
          requirement['tutorials'].each do |tutorial|
            # For each listed tutorial check that a directory with that name exists
            pn = Pathname.new("topics/#{requirement['topic_name']}/tutorials/#{tutorial}")

            if !pn.directory?
              errs.push("Internal requirement to topics/#{requirement['topic_name']}/tutorials/#{tutorial} " \
                        'does not exist')
            end
          end
        end
      when 'none'
        errs.push(*validate_non_empty_key_value(requirement, 'title'))

        requirement.each_key do |x|
          errs.push("Unknown key #{x}") if !%w[title type].include?(x)
        end
      else
        errs.push("Unknown requirement type #{requirement['type']}")
      end
    end

    errs
  end

  def self.lint_faq_file(fn)
    errs = []

    begin
      data = YAML.load_file(fn)
    rescue StandardError => e
      return ["YAML error, failed to parse #{fn}, #{e}"]
    end

    # Check this is something we actually want to process
    if !data.is_a?(Hash)
      puts "Skipping #{fn}"
      return nil
    end
    errs.push(*validate_document(data, @faq_validator))
    errs
  end

  def self.lint_file(fn)
    # Any error messages
    errs = []

    begin
      data = YAML.load_file(fn)
    rescue StandardError => e
      return ["YAML error, failed to parse #{fn}, #{e}"]
      return nil
    end

    # Check this is something we actually want to process
    if !data.is_a?(Hash)
      puts "Skipping #{fn}"
      return nil
    end

    if !fn.include?('metadata.yaml')
      # Load topic metadata for this file
      topic = fn.split('/')[2]
      topic_metadata = YAML.load_file("topics/#{topic}/metadata.yaml")

      # Load subtopic titles
      if data.key?('subtopic')
        subtopic_ids = []
        topic_metadata['subtopics'].each do |x|
          subtopic_ids.push(x['id'])
        end

        @TUTORIAL_SCHEMA['mapping']['subtopic']['enum'] = subtopic_ids
        @SLIDES_SCHEMA['mapping']['subtopic']['enum'] = subtopic_ids
        @tutorial_validator = Kwalify::Validator.new(@TUTORIAL_SCHEMA)
        @slides_validator = Kwalify::Validator.new(@SLIDES_SCHEMA)
      end
    end

    # Generic error handling:
    ## Check requirements
    errs.push(*validate_requirements(data['requirements'])) if data.key?('requirements')

    ## Check follow ups
    errs.push(*validate_requirements(data['follow_up_training'])) if data.key?('follow_up_training')

    # Custom error handling:
    if tutorial?(fn)
      errs.push(*validate_document(data, @tutorial_validator))
    elsif fn.include?('metadata.yaml')
      errs.push(*validate_document(data, @topic_validator))
    elsif slide?(fn)
      errs.push(*validate_document(data, @slides_validator))
    end

    # Check contributors OR contributions
    if (slide?(fn) || tutorial?(fn)) && !(data.key?('contributors') || data.key?('contributions'))
      errs.push('Document lacks EITHER contributors OR contributions key')
    end

    # If we had no errors, validated successfully
    errs
  end

  def self.lint_quiz_file(fn)
    errs = []

    begin
      data = YAML.load_file(fn)
    rescue StandardError => e
      return ["YAML error, failed to parse #{fn}, #{e}"]
      return nil
    end

    if !data.is_a?(Hash)
      warn "Skipping #{fn}"
      return nil
    end

    data['questions'].select { |q| q.key? 'correct' }.each do |q|
      if q['correct'].is_a?(Array)
        if q['type'] != 'choose-many'
          errs.push("There are multiple answers for this question, but it is not a choose-many #{q['title']}")
        end

        q['correct'].each do |c|
          errs.push("Answer #{c} not included in options for question #{q['title']}") if !q['answers'].include?(c)
        end
      else
        if q['type'] != 'choose-1'
          errs.push("There is only a single textual answer, it must be a list for a choose-many question #{q['title']}")
        end

        if !q['answers'].include?(q['correct'])
          errs.push("Answer #{q['correct']} not included in options for question #{q['title']}")
        end
      end
    end

    errs.push(*validate_document(data, @quiz_validator))
    errs
  end

  def self.run
    errors = []
    Find.find('./topics') do |path|
      if FileTest.directory?(path)
        next unless File.basename(path).start_with?('.')

        Find.prune       # Don't look any further into this directory.

      else
        last_component = path.split('/')[-1]
        if last_component =~ /slides.*html$/ || last_component =~ /^tutorial_?(ES|FR|).md/ || last_component =~ /metadata.ya?ml/
          errs = lint_file(path)
          errors += [[path, errs]] if !errs.nil? && errs.length.positive?
        end
      end
    end

    Dir.glob('**/faqs/**/*.md') do |path|
      if FileTest.directory?(path)
        next unless File.basename(path).start_with?('.')

        Find.prune       # Don't look any further into this directory.

      else
        last_component = path.split('/')[-1]
        next if path =~ /aaaa_dontquestionthislinkitisthegluethatholdstogetherthegalaxy.md/

        if last_component =~ (/.*.md/) && !(last_component =~ /index.md$/ || last_component =~ /README.md/)
          errs = lint_faq_file(path)
          if !errs.nil? && errs.length.positive?
            puts "#{path}: #{errs}"
            errors += [[path, errs]]
          end
        end
      end
    end

    Dir.glob('./topics/**/quiz/*') do |path|
      if FileTest.directory?(path)
        next unless File.basename(path).start_with?('.')

        Find.prune       # Don't look any further into this directory.

      else
        last_component = path.split('/')[-1]
        errors += lint_quiz_file(path) if last_component =~ /ya?ml/
      end
    end

    errors
  end
end

if $PROGRAM_NAME == __FILE__
  ec = 0
  errors = SchemaValidator.run
  if errors.length.positive?
    ec = 1
    errors.each do |path, errs|
      # Otherwise, print errors and exit non-zero
      puts "\e[48;5;09m#{path} has errors\e[m"
      errs.each { |x| puts "  #{x}" }
    end
  end
  exit ec
end
