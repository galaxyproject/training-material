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

if __FILE__ == $PROGRAM_NAME
  require 'test/unit'
    # Testing for the class
  class IntersectionTest < Test::Unit::TestCase
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
