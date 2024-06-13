#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'net/http'
require 'csv'
require 'date'
require 'yaml'

# Fetch data from a google sheet
url = 'https://docs.google.com/spreadsheets/d/1iXjLlMEH5QMAMyUMHi1c_Lb7OiJhL_9hgJrtAsBoZ-Y/export?format=tsv'
data = `curl -sL "#{url}"`

data = CSV.parse(data, col_sep: "\t", headers: true, quote_char: '|')
count = 0

# define some columns
col_material = 3
col_length = 5
col_speakers = 6
col_galaxyversion = 10

## recordings metadata definition on tutorials/slides
#
# recordings:
# - speakers:
#    - shiltemann
#    captioners:
#    - hexylena
#    - bebatut
#    date: '2020-06-12'
#    galaxy_version: '20.05'
#    length: 51M
#    youtube_id: "oAVjF_7ensg"

data.each do |row|
  # Parse
  # 29/01/2024 14:04:47
  submission_date = DateTime.strptime(row['Timestamp'], '%d/%m/%Y %H:%M:%S')

  # extract metadata from Google form
  length = row[col_length]
  galaxy_version= row[col_galaxyversion]
  speakers = row[col_speakers]
  date = submission_date.strftime('%Y-%m-%d')

  material_file = row[col_material].gsub("tutorial.html","tutorial.md").gsub("https://training.galaxyproject.org/","").gsub("training-material/","")

  bot_timestamp = submission_date.to_time.to_i
  recording_metadata = {"youtube_id" => "TODO",
                        "length" => length,
                        "galaxy_version" => galaxy_version,
                        "date" => "'#{date}'",
                        "speakers" => "[#{speakers}]",
                        "captioners" => "[#{speakers}]",
                        "bot-timestamp" => bot_timestamp }

  # append metadata into GTN material
  material_metadata = YAML.load_file(material_file)
  puts recording_metadata.to_yaml

  if material_metadata["recordings"]
    # check the "bot_timestamp"
    exists = false
    for rec in material_metadata["recordings"]
      puts rec["bot-timestamp"]
      if rec["bot-timestamp"].to_s == bot_timestamp.to_s
        exists = true
      end
    end

    if exists
      puts "recording already added, skipping"
    else
      puts "new recording, adding"
      material_metadata["recordings"].push(recording_metadata)
    end
  else
    puts "first recording, adding"
    material_metadata["recordings"] = [recording_metadata]
  end

  #pp material_metadata

  # write to file
  material_original = File.open(material_file,"r").read.split("---\n",3)

  outfile = File.open(material_file,"w")
  outfile.write("#{material_metadata.to_yaml}\n\n---\n\n#{material_original[2]}")

end
