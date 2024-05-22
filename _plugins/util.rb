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
  if s.year == e.year
    if s.month == e.month
      if s.day == e.day
        "#{s.strftime('%B')} #{s.day}, #{s.year}"
      else
        "#{s.strftime('%B')} #{s.day}-#{e.day}, #{s.year}"
      end
    else
      "#{s.strftime('%B')} #{s.day}-#{e.strftime('%B')} #{e.day}, #{s.year}"
    end
  else
    "#{s.strftime('%B')} #{s.day}, #{s.year}-#{e.strftime('%B')} #{e.day}, #{e.year}"
  end
end
