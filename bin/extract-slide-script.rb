#!/usr/bin/env ruby
require 'yaml'
require 'json'

fn = ARGV[0]
metadata = YAML.load_file(fn)
WORD_MAP = {}
YAML.load_file(ARGV[1]).each_pair do |k,v|
 WORD_MAP.merge!({k.downcase => v})
end

out_script = File.new(ARGV[2], 'w')

# Do we have these slides? Yes or no.
has_questions = metadata.fetch('questions', []).length > 0
has_objectives = metadata.fetch('objectives', []).length > 0

# Parse the material for the slide notes
file = File.open(fn)
lines = file.readlines.map(&:chomp)

# The structure will be
# ---
# meta
# ---
#
# contents

# +1 because we skipped the 0th entry, +1 again to not include the `---`
end_meta = lines[1..-1].index("---") + 2

# Strip off the metadata
contents = lines[end_meta..-1]

# This will be our final script
blocks = []
if has_questions
  blocks.push(metadata['questions'])
end
if has_objectives
  blocks.push(metadata['objectives'])
end

# Accumulate portions between ??? and ---
current_block = []
in_notes = false
contents.each{ |x|
  # Check whether we're in the notes or out of them.
  if x == "???" then
    in_notes = true
  elsif x == "---" then
    if in_notes then
      blocks.push(current_block)
      current_block = []
    end

    in_notes = false
  end

  if in_notes then
    current_block.push(x)
  end
}
blocks.push(current_block)

def translate(word)
  m = /([^A-Za-z0-9]*)([A-Za-z0-9]+)([^A-Za-z0-9]*)(.*)/.match(word)


  if m[2] then
    fixed = WORD_MAP.fetch(m[2].downcase, m[2])
  else
    fixed = m[2]
  end

  #puts "#{m} ⇒ #{m[1] + fixed + m[3]}"
  return m[1] + fixed + m[3] + m[4]
end

# For each block, cleanup first.
blocks = blocks.map{ |block|
  # Remove the - prefix from each line
  script_lines = block.map{ |x| x.delete_prefix("- ") }
  # Remove the leading ???
  script_lines = script_lines[1..-1]
  # Remove blank entries
  script_lines = script_lines.filter{ |x| x.length != 0 }
  # Translate specific words as needed
  script_lines = script_lines.map{ |line|
    line.delete_prefix("- ")
    # If they don't end with punctuation, fix it.
    if ! (line.end_with?('.') or line.end_with?('?') or line.end_with?('!'))
      line += '.'
    end

    line
  }
  script_lines
}

#out_subs.write(blocks.map{ |line| line.join(" ") }.join("\n"))

blocks2 = blocks.map { |block|
  s = block.map{ |line|
    line = line.split(' ').map{ |w|
      translate(w)
    }.join(' ')
    line
  }
  s
}

final = blocks.zip(blocks2)
dump = JSON.pretty_generate(final)
out_script.write(dump)
