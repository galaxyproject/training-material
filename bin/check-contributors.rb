#!/usr/bin/env ruby
# frozen_string_literal: true

# Checks the header of a tutorial or slides file fo the current contributor list
# Then compares against the people who have touched that file in git using the git log
#
# Each tutorial or slides file needs to be checked individually.
#
# If there are unknown contributors, create a .mailmap file like:
# <gh-username> <user@earlham.ac.uk>
# <gh-username> <user@gmail.com>
#
# This is the format of git's mail-map
# https://www.git-scm.com/docs/git-check-mailmap
# We probably should not commit this file, some people don't want their old emails
# listed so publicly.

require 'yaml'

if ARGV.size != 1
  puts 'Please run with ./bin/check-contributors path/to/tutorial.md/or/slides.html'
  exit
end

fn = ARGV[0]

# Any error messages
# errs = []
data = YAML.load_file(fn)
current_contributors = data['contributors']

# Full Contributors Data
CONTRIBUTORS = YAML.load_file('CONTRIBUTORS.yaml')
contributor_emails = CONTRIBUTORS.map do |k, v|
  [v['email'], k] if v && v.key?('email')
end.compact.to_h

file_contributors = `git log --use-mailmap --follow --pretty=%aE #{fn}`.lines.sort.uniq

fixed_contribs = file_contributors.map do |email|
  email = email.strip
  if /users.noreply.github.com/.match(email)
    parts = /^(?<_num>[0-9]+\+)?(?<id>.*)@users.noreply.github.com/.match(email)
    # we just want their gh id
    parts[:id]
  elsif contributor_emails.key?(email)
    contributor_emails[email]
  else
    email
  end
end

# known contributors
known = fixed_contribs.reject { |x| /@/.match(x) }
unknown = fixed_contribs.select { |x| /@/.match(x) }

missing = (known - current_contributors).sort.uniq
# These contributors not yet recognised
puts "Missing contributors: #{missing}" if missing.length.positive?

# These emails might map to known users, but we don't know yet.
puts "Unknown emails: #{unknown.sort.uniq}" if unknown.length.positive?
