#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'net/http'
require 'csv'
require 'date'

# Fetch data from a google sheet
url = 'https://docs.google.com/spreadsheets/d/1KJsQpZsUfaZh-a8jLx0gNymaQxvSNForrekEgIzF3Eo/export?format=tsv'
data = `curl -sL "#{url}"`

contributions_keys = %w[Authorship Editing Testing Infrastructure Translation Funding]

data = CSV.parse(data, col_sep: "\t", headers: true, quote_char: '|')
count = 0

data.each do |row|
  # Parse
  # 29/01/2024 14:04:47
  post_date = DateTime.strptime(row['Timestamp'], '%d/%m/%Y %H:%M:%S')
  filename = "news/_posts/#{post_date.strftime('%Y-%m-%d')}-#{row['Title'].downcase.gsub(/[^a-z0-9\s]/i, '').gsub(
    /\s+/, '-'
  )}.md"

  # Skip some testing posts
  if (row['Title'] == 'TESTING') || (row['Title'] == "Wendi's Dog is the Best")
    STDERR.puts "Skipping #{filename} as it is a test post"
    next
  end

  # Don't overwrite existing posts
  if File.exist?(filename)
    STDERR.puts "Skipping #{filename} as it already exists"
    next
  end

  STDERR.puts "Creating #{filename}"
  count += 1

  post_metadata = {
    'title' => row['Title'],
    'layout' => 'news',
    'tags' => row['Tags'].split(',').map(&:strip),
    'from_google_form' => true,
    'contributions' => {}
  }

  contributions_keys.each do |key|
    post_metadata['contributions'][key.downcase] = row[key].split(',').map(&:strip) if row[key]
  end

  if row['Optional Cover Image']
    post_metadata['cover'] = row['Optional Cover Image']
    post_metadata['coveralt'] = row['Cover Image Alternative Text']
  end

  if row['Link to a specific tutorial you wish to push people to view']
    post_metadata['tutorial'] = row['Link to a specific tutorial you wish to push people to view']
  end

  if row['Link for a non-tutorial thing you want to push people to view']
    post_metadata['link'] = row['Link for a non-tutorial thing you want to push people to view']
  end

  # Serialise to a file

  File.open(filename, 'w') do |file|
    file.puts YAML.dump(post_metadata)
    file.puts "---\n"
    file.puts row['Blog Post'].gsub('  ', "\n\n")
  end
end

puts count
