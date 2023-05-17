#!/usr/bin/env ruby
# frozen_string_literal: true

require 'yaml'
require 'shellwords'
require 'json'
require 'find'
require 'bibtex'
require 'citeproc/ruby'
require 'csl/styles'

fn = ARGV[0]
metadata = YAML.load_file(fn)

topic_fn = "#{fn.split('/').slice(0, 2).join('/')}/metadata.yaml"
topic_metadata = YAML.load_file(topic_fn)

ARI_MAP = File.expand_path(File.join(__dir__, 'ari-map.yml'))
WORD_MAP = {}
YAML.load_file(ARI_MAP).each_pair do |k, v|
  WORD_MAP.merge!({ k.downcase => v })
end

APPROVED_VOICES = {
  'en' => [
    { 'id' => 'Amy', 'lang' => 'en-GB', 'neural' => true },
    { 'id' => 'Aria', 'lang' => 'en-NZ', 'neural' => true },
    { 'id' => 'Brian', 'lang' => 'en-GB', 'neural' => true },
    { 'id' => 'Emma', 'lang' => 'en-GB', 'neural' => true },
    { 'id' => 'Joanna', 'lang' => 'en-US', 'neural' => true },
    { 'id' => 'Joey', 'lang' => 'en-US', 'neural' => true },
    { 'id' => 'Kendra', 'lang' => 'en-US', 'neural' => true },
    { 'id' => 'Matthew', 'lang' => 'en-US', 'neural' => true },
    { 'id' => 'Nicole', 'lang' => 'en-AU', 'neural' => false },
    { 'id' => 'Olivia', 'lang' => 'en-AU', 'neural' => true },
    { 'id' => 'Raveena', 'lang' => 'en-IN', 'neural' => false },
    { 'id' => 'Salli', 'lang' => 'en-US', 'neural' => true },
    { 'id' => 'Ayanda', 'lang' => 'en-ZA', 'neural' => true },
    { 'id' => 'Geraint', 'lang' => 'en-GB-WLS', 'neural' => false }
  ],
  'es' => [
    { 'id' => 'Miguel', 'lang' => 'es-US', 'neural' => false },
    { 'id' => 'Mia', 'lang' => 'es-MX', 'neural' => false },
    { 'id' => 'Enrique', 'lang' => 'es-ES', 'neural' => false },
    { 'id' => 'Conchita', 'lang' => 'es-ES', 'neural' => false },
    { 'id' => 'Lupe', 'lang' => 'es-US', 'neural' => true }
  ]
}.freeze

# This is copied directly from the plugins, TODO: make into a module.
global_bib = BibTeX::Bibliography.new
bib_paths = [Find.find('./topics'), Find.find('./faqs')].lazy.flat_map(&:lazy)
bib_paths.each  do |path|
  if FileTest.directory?(path)
    next unless File.basename(path).start_with?('.')

    Find.prune # Don't look any further into this directory.

  elsif path =~ /bib$/
    BibTeX.open(path).each do |x|
      x = x.convert_latex
      global_bib << x
    end
  end
end
cp = CiteProc::Processor.new format: 'text', locale: 'en'
cp.import global_bib.to_citeproc

# Do we have these slides? Yes or no.
m_qs = metadata.fetch('questions', [])
m_qs = [] if m_qs.nil?
has_questions = m_qs.length.positive?

m_os = metadata.fetch('objectives', [])
m_os = [] if m_os.nil?
has_objectives = m_os.length.positive?

m_kp = metadata.fetch('key_points', [])
m_kp = [] if m_kp.nil?
has_keypoints = m_kp.length.positive?

m_rq = metadata.fetch('requirements', [])
m_rq = [] if m_rq.nil?
t_rq = topic_metadata.fetch('requirements', [])
t_rq = [] if t_rq.nil?
has_requirements = m_rq.length.positive? || t_rq.length.positive?

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
end_meta = lines[1..].index('---') + 2

# Strip off the metadata
contents = lines[end_meta..]

# This will be our final script
blocks = [[metadata['title']]]
if has_requirements
  if m_lang == 'es'
    blocks.push(['Antes de profundizar en el contenido de estas diapositivas, te recomendamos que le des un vistazo a'])
  else
    blocks.push(['Before diving into this slide deck, we recommend you to have a look at the following.'])
  end
end
blocks.push(metadata['questions']) if has_questions
blocks.push(metadata['objectives']) if has_objectives

# Accumulate portions between ??? and ---
current_block = []
in_notes = false
contents.each do |x|
  # Check whether we're in the notes or out of them.
  if x == '???'
    in_notes = true
  elsif ['---', '--'].include?(x)
    if in_notes
      blocks.push(current_block)
      current_block = []
    end

    in_notes = false
  end

  current_block.push(x) if in_notes
end
blocks.push(current_block)
blocks.push(metadata['key_points']) if has_keypoints

if m_lang == 'es'
  blocks.push(['¡Gracias por ver este vídeo!'])
else
  blocks.push(['Thank you for watching!'])
end

# For each block, cleanup first.
blocks = blocks.map do |block|
  # Remove the - prefix from each line
  script_lines = block.map { |x| x.strip.delete_prefix('- ') }
  # Remove the leading ???
  script_lines = script_lines[1..] if script_lines[0] == '???'
  # Remove blank entries
  script_lines = script_lines.reject(&:empty?)
  script_lines = script_lines.map do |line|
    line.delete_prefix('- ')
    line.gsub!(/`/, '"')
    # If they don't end with punctuation, fix it.
    line += '.' if !(line.end_with?('.') || line.end_with?('?') || line.end_with?('!'))

    line
  end
  script_lines = script_lines.map do |line|
    line.gsub!(/{%\s*cite ([^}]*)\s*%}/) do |match|
      # Strip off the {% %} first, whitespace, and then remove cite at the
      # start and restrip again.
      value = match[2..-3].strip[4..].strip
      # Render the citation, the :text format includes ( ) on both sides which
      # we strip off.
      cp.render(:citation, id: value)[1..-2]
    end
    line
  end
  script_lines
end

# out_subs.write(blocks.map{ |line| line.join(" ") }.join("\n"))
res = {}
res['blocks'] = blocks

res['voice'] = if m_voice.nil?
                 if m_lang == 'es'
                   APPROVED_VOICES['es'].sample
                 else
                   APPROVED_VOICES['en'].sample
                 end
               else
                 metadata['voice']
               end

print JSON.pretty_generate(res)
