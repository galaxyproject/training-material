require 'yaml'

# This module specifically avoids having any dependencies. If you have a method
# that needs to be used in many places, even those without initialisation, it
# should go here.

ALLOWED_SHORT_IDS = [
  'ChangeCase',
  'Convert characters1',
  'Count1',
  'Cut1',
  'Extract_features1',
  'export_remote',
  'Filter1',
  'Grep1',
  'Grouping1',
  'Paste1',
  'Remove beginning1',
  'Show beginning1',
  'Summary_Statistics1',
  'addValue',
  'cat1',
  'comp1',
  'gene2exon1',
  'gff2bed1',
  'intermine',
  'join1',
  'param_value_from_file',
  'random_lines1',
  'sort1',
  'csv_to_tabular',
  # 'ucsc_table_direct1', # This does not work, surprisingly.
  'upload1',
  'wc_gnu',
  'wig_to_bigWig'
].freeze

ALLOWED_LOWER_SHORT_IDS = ALLOWED_SHORT_IDS.map(&:downcase)

def safe_load_yaml(file)
  YAML.load_file(file)
rescue StandardError
  YAML.load_file(file, permitted_classes: [Date])
end

def collapse_event_date_pretty(event)
  s = event['date_start']
  e = if event['date_end'].nil?
        s
      else
        event['date_end']
      end
  # want dates like "Mar 22-25, 2024" or "Mar 22-May 1, 2024"
  dash = ' – ' # thin space, en dash, thin space
  if s.year == e.year
    if s.month == e.month
      if s.day == e.day
        "#{s.strftime('%B')} #{s.day}, #{s.year}"
      else
        "#{s.strftime('%B')} #{s.day}#{dash}#{e.day}, #{s.year}"
      end
    else
      "#{s.strftime('%B')} #{s.day}#{dash}#{e.strftime('%B')} #{e.day}, #{s.year}"
    end
  else
    "#{s.strftime('%B')} #{s.day}, #{s.year}#{dash}#{e.strftime('%B')} #{e.day}, #{e.year}"
  end
end

def safe_site_config(site, key, default)
  if !site.config.nil? && site.config.key?(key)
    site.config[key]
  else
    default
  end
end


def url_prefix(site)
  if !site.config.nil? && site.config.key?('url')
    "#{site.config['url']}#{site.config['baseurl']}"
  else
    'http://localhost:4000/training-material/'
  end
end

def markdownify(site, text)
  site.find_converter_instance(
    Jekyll::Converters::Markdown
  ).convert(text.to_s)
end

def unsafe_slugify(text)
  text.gsub(%r{["'\\/;:,.!@#$%^&*()]}, '').gsub(/\s/, '-').gsub(/-+/, '-')
end

def fix_version(version)
  version
  .gsub('_beta+galaxy', '+galaxy')
  .gsub(/^([0-9]+)_([0-9]+)_([0-9]+)\+galaxy(.+)/, '\1.\2.\3galaxy\4')
  .gsub(/^([0-9]+)\+galaxy(.+)/, '\1.0.0galaxy\2')
  .gsub(/^([0-9.]+)_([0-9]+)/, '\1galaxy\2')
  .gsub(/_rc(.*)galaxy/, 'rc\1galaxy')
  .gsub('+', '')
  .gsub(/^v/, '')
end

def acceptable_tool?(tool_id)
  if ! tool_id.is_a?(String)
    return false
  end

  # Public TS links are fine
  if tool_id.start_with?('toolshed.g2.bx.psu.edu')
    return true
  end

  # These are always allowed (mostly built-ins)
  if ALLOWED_LOWER_SHORT_IDS.include?(tool_id.downcase) || tool_id =~ /^__.*__$/
    return true
  end

  if tool_id.start_with?('interactive_tool_')
    return true
  end

  if tool_id.start_with?('CONVERTER_')
    return true
  end

  # Templated tool IDs are hopefully fine!
  if tool_id.start_with?('{{')
    return true
  end

  # Just the tutorial
  if tool_id.start_with?('Toolshed ID')
    return true
  end

  # Unacceptable
  return false
end


def tool_id_extractor(wf, path: [])
  res = []
  wf['steps'].each do |step_id, step|
    if step.key?('subworkflow')
      res += tool_id_extractor(step['subworkflow'], path: path + [step_id])
    elsif step.key?('tool_id') && ! step['tool_id'].nil?
      res.push(["#{path.join('/')}/#{step_id}", step['tool_id']])
    end
  end
  res
end

def objectify(attrs, url, path)
  obj = attrs.clone
  obj['__path'] = path
  obj['__url'] = url

  def obj.data
    self
  end

  def obj.path
    self['__path']
  end

  def obj.url
    self['__url']
  end

  def obj.content
    self.fetch('content', 'NO CONTENT AVAILABLE')
  end

  def obj.title
    self['title']
  end

  obj
end

if __FILE__ == $PROGRAM_NAME
  require 'test/unit'
    # Testing for the class
  class Gtn::Test::IntersectionTest < Test::Unit::TestCase
    def test_bad_versions
      # toolshed.g2.bx.psu.edu/repos/wolma/mimodd_main/mimodd_info/0.1.8_1
      assert_equal(fix_version("0.1.8_1"), "0.1.8galaxy1")

      # toolshed.g2.bx.psu.edu/repos/iuc/snap_training/snap_training/2013_11_29+galaxy1
      assert_equal(fix_version("2013_11_29+galaxy1"), "2013.11.29galaxy1")

      # toolshed.g2.bx.psu.edu/repos/devteam/vcffilter/vcffilter2/1.0.0_rc1+galaxy3
      assert_equal(fix_version("1.0.0_rc1+galaxy3"), "1.0.0rc1galaxy3")

      #
      assert_equal(fix_version("3+galaxy0"), "3.0.0galaxy0")
    end
  end
end
