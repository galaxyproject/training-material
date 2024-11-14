require 'yaml'

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
