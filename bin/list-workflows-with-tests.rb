#!/usr/bin/env ruby
# frozen_string_literal: true

require 'json'
require 'yaml'

Dir.glob('./topics/**/*.ga') do |path|
  folder = File.dirname(path)
  basename = File.basename(path).gsub(/.ga$/, '')
  possible_tests = Dir.glob("#{folder}/#{basename}*ym*")
  possible_tests = possible_tests.grep(/#{basename}[_-]tests?.ya?ml/)

  possible_tests.each do |possib|
    if !possib.match(/-test.yml/)
      puts "Renaming #{possib}"
      # Rename the file to have the correct extension
      File.rename(possib, possib.gsub(/[_-]tests?.ya?ml$/, '-test.yml'))
    end
  end

  puts path if !possible_tests.empty?
end
