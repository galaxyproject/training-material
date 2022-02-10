#!/usr/bin/env ruby
require 'yaml'
require 'shellwords'
require 'json'

fn = ARGV[0]
metadata = YAML.load_file(fn)

topic_fn = fn.split('/').slice(0, 2).join('/') + '/metadata.yaml'
topic_metadata = YAML.load_file(topic_fn)

ARI_MAP = File.expand_path(File.join(__dir__, 'ari-map.yml'))
WORD_MAP = {}
YAML.load_file(ARI_MAP).each_pair do |k,v|
 WORD_MAP.merge!({k.downcase => v})
end

APPROVED_VOICES = {
  "en" => [
    {"id" =>"Amy"     , "lang" => "en-GB" , "neural" => true},
    {"id" =>"Aria"    , "lang" => "en-NZ" , "neural" => true},
    {"id" =>"Brian"   , "lang" => "en-GB" , "neural" => true},
    {"id" =>"Emma"    , "lang" => "en-GB" , "neural" => true},
    {"id" =>"Joanna"  , "lang" => "en-US" , "neural" => true},
    {"id" =>"Joey"    , "lang" => "en-US" , "neural" => true},
    {"id" =>"Kendra"  , "lang" => "en-US" , "neural" => true},
    {"id" =>"Matthew" , "lang" => "en-US" , "neural" => true},
    {"id" =>"Nicole"  , "lang" => "en-AU" , "neural" => false},
    {"id" =>"Olivia"  , "lang" => "en-AU" , "neural" => true},
    {"id" =>"Raveena" , "lang" => "en-IN" , "neural" => false},
    {"id" =>"Salli"   , "lang" => "en-US" , "neural" => true}
  ],
  "es" => [
    { "id" => "Miguel"   , "lang" => "es-US" , "neural" => false },
    { "id" => "Mia"      , "lang" => "es-MX" , "neural" => false },
    { "id" => "Enrique"  , "lang" => "es-ES" , "neural" => false },
    { "id" => "Conchita" , "lang" => "es-ES" , "neural" => false },
    { "id" => "Lupe"     , "lang" => "es-US" , "neural" => true }
  ]
}

# Do we have these slides? Yes or no.
m_qs = metadata.fetch('questions', [])
m_qs = [] if m_qs.nil?
has_questions = m_qs.length > 0

m_os = metadata.fetch('objectives', [])
m_os = [] if m_os.nil?
has_objectives = m_os.length > 0

m_kp = metadata.fetch('key_points', [])
m_kp = [] if m_kp.nil?
has_keypoints = m_kp.length > 0

m_rq = metadata.fetch('requirements', [])
m_rq = [] if m_rq.nil?
t_rq = topic_metadata.fetch('requirements', [])
t_rq = [] if t_rq.nil?
has_requirements = m_rq.length > 0 || t_rq.length > 0

m_lang = metadata.fetch('lang', 'en')
m_voice = metadata.fetch('voice', nil)

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
blocks = [[metadata['title']]]
if has_requirements
  if m_lang == "en" then
    blocks.push(['Before diving into this slide deck, we recommend you to have a look at the following.'])
  elsif m_lang == "es" then
    blocks.push(["Antes de profundizar en el contenido de estas diapositivas, te recomendamos que le des un vistazo a"])
  else
    blocks.push(['Before diving into this slide deck, we recommend you to have a look at the following.'])
  end
end
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
  elsif x == "---" or x == "--" then
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
if has_keypoints
  blocks.push(metadata['key_points'])
end

if m_lang == "en" then
  blocks.push(["Thank you for watching!"])
elsif m_lang == "es" then
  blocks.push(["¡Gracias por ver este vídeo!"])
else
  blocks.push(["Thank you for watching!"])
end

# For each block, cleanup first.
blocks = blocks.map{ |block|
  # Remove the - prefix from each line
  script_lines = block.map{ |x| x.strip.delete_prefix("- ") }
  # Remove the leading ???
  if script_lines[0] == '???'
    script_lines = script_lines[1..-1]
  end
  # Remove blank entries
  script_lines = script_lines.select{ |x| x.length != 0 }
  script_lines = script_lines.map{ |line|
    line.delete_prefix("- ")
    line.gsub!(/`/, '"')
    # If they don't end with punctuation, fix it.
    if ! (line.end_with?('.') or line.end_with?('?') or line.end_with?('!'))
      line += '.'
    end

    line
  }
  script_lines
}

#out_subs.write(blocks.map{ |line| line.join(" ") }.join("\n"))
res = Hash.new
res["blocks"] = blocks

if m_voice.nil? then
  if m_lang == "en" then
    res["voice"] = APPROVED_VOICES['en'].sample
  elsif m_lang == "es" then
    res["voice"] = APPROVED_VOICES['es'].sample
  else
    res["voice"] = APPROVED_VOICES['en'].sample
  end
else
  res["voice"] = metadata['voice']
end

print JSON.pretty_generate(res)
