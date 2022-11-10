#!/usr/bin/env ruby
require 'yaml'
require 'find'
require 'pathname'
require 'kwalify'
require './bin/gtn.rb'


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
  $topic_validator = Kwalify::Validator.new(@TOPIC_SCHEMA)
  $tutorial_validator = Kwalify::Validator.new(@TUTORIAL_SCHEMA)
  $slides_validator = Kwalify::Validator.new(@SLIDES_SCHEMA)
  $faq_validator = Kwalify::Validator.new(@FAQ_SCHEMA)
  $quiz_validator = Kwalify::Validator.new(@QUIZ_SCHEMA)
  $requirement_external_validator = Kwalify::Validator.new(@requirement_external_schema)
  $requirement_internal_validator = Kwalify::Validator.new(@requirement_internal_schema)

  def self.validate_document(document, validator)
    errors = validator.validate(document)
    if errors && !errors.empty?
      return errors
    end
    return []
  end

  def self.validate_non_empty_key_value(map, key)
      if map.key?(key) then
        if map[key].length == 0 then
          return ["Empty #{key} for requirement"]
        end
      else
        return ["Missing #{key} for requirement"]
      end
      return []
  end

  def self.is_tutorial(fn)
    fn.include?('tutorial.md') || fn =~ /tutorial_[A-Z]{2,}.md/
  end

  def self.is_slide(fn)
    fn.include?('slides.html') || fn =~ /slides_[A-Z]{2,}.html/
  end

  def self.validate_requirements(requirements)
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

  def self.lint_faq_file(fn)
    errs = []

    begin
      data = YAML.load_file(fn)
    rescue
      puts "Skipping #{fn}"
      return nil
    end

    # Check this is something we actually want to process
    if ! data.is_a?(Hash) then
      puts "Skipping #{fn}"
      return nil
    end
    errs.push(*validate_document(data, $faq_validator))
    return errs
  end

  def self.lint_file(fn)
    # Any error messages
    errs = []

    begin
      data = YAML.load_file(fn)
    rescue
      puts "Skipping #{fn}"
      return nil
    end

    # Check this is something we actually want to process
    if ! data.is_a?(Hash) then
      puts "Skipping #{fn}"
      return nil
    end

    if not fn.include?('metadata.yaml') then
      # Load topic metadata for this file
      topic = fn.split('/')[2]
      topic_metadata = YAML.load_file("topics/#{topic}/metadata.yaml")

      # Load subtopic titles
      if data.key?('subtopic') then
        subtopic_ids = []
        topic_metadata['subtopics'].each{ |x|
          subtopic_ids.push(x['id'])
        }

        @TUTORIAL_SCHEMA['mapping']['subtopic']['enum'] = subtopic_ids
        @SLIDES_SCHEMA['mapping']['subtopic']['enum'] = subtopic_ids
      $tutorial_validator = Kwalify::Validator.new(@TUTORIAL_SCHEMA)
      $slides_validator = Kwalify::Validator.new(@SLIDES_SCHEMA)
      end
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
    if is_tutorial(fn) then
      errs.push(*validate_document(data, $tutorial_validator))
    elsif fn.include?('metadata.yaml') then
      errs.push(*validate_document(data, $topic_validator))
    elsif is_slide(fn) then
      errs.push(*validate_document(data, $slides_validator))
    else
      #errs.push("No validation available for this type of file")
    end

    # Check contributors OR contributions
    if is_slide(fn) || is_tutorial(fn) then
      if not (data.has_key?('contributors') or data.has_key?('contributions'))
        errs.push("Document lacks EITHER contributors OR contributions key")
      end
    end

    # If we had no errors, validated successfully
    return errs
  end

  def self.lint_quiz_file(fn)
    errs = []

    begin
      data = YAML.load_file(fn)
    rescue
      $stderr.puts "Skipping #{fn}"
      return nil
    end

    if ! data.is_a?(Hash) then
      $stderr.puts "Skipping #{fn}"
      return nil
    end

    data['questions'].select{|q| q.has_key? "correct"}.each{|q|
      if q["correct"].is_a?(Array) then
        if q["type"] != "choose-many" then
          errs.push("There are multiple answers for this question, but it is not a choose-many #{q['title']}")
        end

        q["correct"].each{|c|
          if ! q["answers"].include?(c)
            errs.push("Answer #{c} not included in options for question #{q['title']}")
          end
        }
      else
        if q["type"] != "choose-1" then
          errs.push("There is only a single textual answer, it must be a list for a choose-many question #{q['title']}")
        end

        if ! q["answers"].include?(q["correct"])
          errs.push("Answer #{q["correct"]} not included in options for question #{q['title']}")
        end
      end
    }

    errs.push(*validate_document(data, $quiz_validator))
    return errs
  end

  def self.run
    errors = []
    Find.find('./topics') do |path|
      if FileTest.directory?(path)
        if File.basename(path).start_with?('.')
          Find.prune       # Don't look any further into this directory.
        else
          next
        end
      else
        last_component = path.split('/')[-1]
        if last_component =~ /slides.*html$/ || last_component =~ /tutorial.*md/  || last_component =~ /metadata.ya?ml/  then
          errs = lint_file(path)
          if ! errs.nil? and errs.length > 0
            errors += [ [path, errs] ]
          end
        end
      end
    end


    Dir.glob("**/faqs/**/*.md") do |path|
      if FileTest.directory?(path)
        if File.basename(path).start_with?('.')
          Find.prune       # Don't look any further into this directory.
        else
          next
        end
      else
        last_component = path.split('/')[-1]
        if path =~ /aaaa_dontquestionthislinkitisthegluethatholdstogetherthegalaxy.md/
          next
        end
        if last_component =~ /.*.md/ then
          unless last_component =~ /index.md$/ || last_component =~ /README.md/  then
            errs = lint_faq_file(path)
            if ! errs.nil? and errs.length > 0
              puts "#{path}: #{errs}"
              errors += [ [path, errs] ]
            end
          end
        end
      end
    end

    Dir.glob("./topics/**/quiz/*") do |path|
      if FileTest.directory?(path)
        if File.basename(path).start_with?('.')
          Find.prune       # Don't look any further into this directory.
        else
          next
        end
      else
        last_component = path.split('/')[-1]
        if last_component =~ /ya?ml/  then
          errors += lint_quiz_file(path)
        end
      end
    end

    errors
  end
end

if $0 == __FILE__
  ec = 0
  errors = SchemaValidator::run()
  if errors.length > 0
    ec = 1
    errors.each{|path, errs|
      # Otherwise, print errors and exit non-zero
      puts "\e[48;5;09m#{path} has errors\e[m"
      errs.each {|x| puts "  #{x}" }
    }
  end
  exit ec
end
