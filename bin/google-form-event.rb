#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'net/http'
require 'csv'
require 'date'
require './_plugins/util'

# Fetch data from a google sheet
url = 'https://docs.google.com/spreadsheets/d/1LTShstcORXf_zB06naOfPzisAoOcik-xJvgK5LQCGP8/export?format=tsv'
data = `curl -sL "#{url}"`

data = CSV.parse(data, col_sep: "\t", headers: true, quote_char: '|')
count = 0
pr_message = ""

data.each do |row|
  # Parse
  # 29/01/2024 14:04:47
  event_date = DateTime.strptime(row['Start date'], '%d/%m/%Y')
  if row['End date']
    event_date_end = DateTime.strptime(row['End date'], '%d/%m/%Y')
  end
  post_date = DateTime.strptime(row['Timestamp'], '%d/%m/%Y %H:%M:%S')

  filename = "events/#{event_date.strftime('%Y-%m-%d')}-#{row['Title of your Event'].downcase.gsub(/[^a-z0-9\s-]/i, '').gsub(/\s+/, ' ').gsub(/ /, '-')}.md"

  # Skip some testing posts
  if (row['Title of your Event'] == 'TESTING')
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
    'layout' => 'event-external',
    'google_form_id' => post_date.to_time.to_i,
    'title' => row['Title of your Event'],
    'description' => row['Description of your event'],
    'external' => row['Link to event page'],
    'contributions' => {
       'organisers' => row['Organizers already in the GTN CONTRIBUTORS file'].split(',').map(&:strip)
    },
    'location' => {
        'name' => row['Location of the event']
    },
    'date_start' => event_date.to_date,
  }

  if row['End date']
    post_metadata['date_end'] = event_date_end.to_date
  end

  # Serialise to a file
  File.open(filename, 'w') do |file|
    file.puts YAML.dump(post_metadata)
    file.puts "---\n"

    if row['Comments'] && row['Comments'].length.positive?
      file.puts "<!-- #{row['Comments']} -->"
    end
  end


  if row['Organizers not yet in the GTN CONTRIBUTORS file']
    pr_message += "<br>TODO: add the following contributors to the GTN: <br> #{row['Organizers not yet in the GTN CONTRIBUTORS file']}"
  end
  if row['Anything else we should know?']
    pr_message += "<br><br>Remarks from submitter:<br> #{row['Anything else we should know?']}"
  end
end


puts "new_ids=#{count}"
puts "pr_message=#{pr_message}"
