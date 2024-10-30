#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'net/http'
require 'csv'
require 'date'
require './_plugins/util'

# Fetch data from a google sheet
url = 'https://docs.google.com/spreadsheets/d/1RFF3G9_bP8EpfACBk8lnF-Ib43ZGGAMm3ewPwW7eFT0/export?format=tsv'
data = `curl -sL "#{url}"`

data = CSV.parse(data, col_sep: "\t", headers: true)
count = 0

data.each do |row|
  # Parse
  # 29/01/2024 14:04:47
  post_date = DateTime.strptime(row['Timestamp'], '%d/%m/%Y %H:%M:%S')
  filename = "faqs/#{row['This FAQ Concerns'].downcase}/#{row['Title'].downcase.gsub(/[^a-z0-9\s-]/i, '').gsub(/\s+/, ' ').gsub(/ /, '-')}.md"

  # Skip some testing posts
  if (row['Title'] == 'TESTING')
    STDERR.puts "Skipping #{filename} as it is a test post"
    next
  end

  # Don't overwrite existing posts
  if File.exist?(filename)
    other_file = safe_load_yaml(filename)
    if other_file['google_form_id'] == post_date.to_time.to_i
      STDERR.puts "Skipping #{filename} as it already exists"
      next
    else
      filename = filename.gsub('.md', "-#{post_date.to_time.to_i}.md")
    end
  end

  STDERR.puts "Creating #{filename}"
  count += 1

  post_metadata = {
    'title' => row['Title'],
    'layout' => 'faq',
    'area' => row['Area'],
    'box_type' => 'tip',
    'google_form_id' => post_date.to_time.to_i,
    'contributors' => [row['Your GitHub username']]
  }

  # Serialise to a file
  File.open(filename, 'w') do |file|
    file.puts YAML.dump(post_metadata)
    file.puts "---\n"
    file.puts row['FAQ Text'].gsub('  ', "\n\n")

    if row['Comments'] && row['Comments'].length.positive?
      file.puts "<!-- #{row['Comments']} -->"
    end
  end
end

puts count
