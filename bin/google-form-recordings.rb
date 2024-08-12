#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'net/http'
require 'csv'
require 'date'
require './_plugins/util'

# Fetch data from a google sheet
url = 'https://docs.google.com/spreadsheets/d/1iXjLlMEH5QMAMyUMHi1c_Lb7OiJhL_9hgJrtAsBoZ-Y/export?format=tsv'
data = `curl -sL "#{url}"`
new_recordings = false

data = CSV.parse(data, col_sep: "\t", headers: true, quote_char: '|')
count = 0

# define some columns
col_material = 3
col_length = 5
col_speakers = 6
col_galaxyversion = 10
col_prmade = 12

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
  speakers = row[col_speakers].split(",")
  date = submission_date.strftime('%Y-%m-%d')

  if row[col_material] == 'TESTING' or row[col_prmade] == 'yes'
    STDERR.puts "Skipping recording as it is a test or a PR was already openened"
    next
  end

  material_file = row[col_material].gsub("tutorial.html","tutorial.md").gsub("https://training.galaxyproject.org/","").gsub("training-material/","").gsub("training-material//","")


  bot_timestamp = submission_date.to_time.to_i
  recording_metadata = {"youtube_id" => "TODO",
                        "length" => length,
                        "galaxy_version" => galaxy_version,
                        "date" => date,
                        "speakers" => speakers,
                        "captioners" => speakers.dup,
                        "bot-timestamp" => bot_timestamp }

  # append metadata into GTN material
  material_metadata = safe_load_yaml(material_file)

  if material_metadata["recordings"]
    # check the "bot_timestamp"
    exists = false
    for rec in material_metadata["recordings"]
      if rec["bot-timestamp"].to_s == bot_timestamp.to_s
        exists = true
      end
    end

    if !exists
      material_metadata["recordings"].push(recording_metadata)
      new_recordings = true
    end
  else
    material_metadata["recordings"] = [recording_metadata]
    new_recordings = true
  end

  #pp material_metadata

  # write to file
  material_original = File.open(material_file,"r").read.split("---\n",3)

  outfile = File.open(material_file,"w")
  outfile.write("#{material_metadata.to_yaml}\n\n---\n\n#{material_original[2]}")

end

STDERR.puts "new recordings: #{new_recordings}"
puts new_recordings
