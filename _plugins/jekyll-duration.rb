module Jekyll
  module DurationFilter
    def duration_to_human(duration)
      # Match the different parts of the string, must match entire string or it will fail.
      match = /^(?:([0-9]*)[Hh])*(?:([0-9]*)[Mm])*(?:([0-9.]*)[Ss])*$/.match(duration)

      # If it doesn't match, pass through unedited so we don't cause unexpected issues.
      if match.nil? then
        puts "Could not parse time: #{duration}"
        return duration
      end

      # Otherwise append english terms for the various parts
      duration_parts = []

      hour = "hour"
      hours = "hours"
      minutes = "minutes"
      if @context.registers[:page]&.key?('lang')
          lang = @context.registers[:page]['lang']
          hour = @context.registers[:site].data["lang"][lang]["hour"]
          hours = @context.registers[:site].data["lang"][lang]["hours"]
          minutes = @context.registers[:site].data["lang"][lang]["minutes"]
      end

      # Hours
      if ! match[1].nil? then
        if match[1] == '1' then
          duration_parts.push("#{match[1]} "+hour)
        else
          duration_parts.push("#{match[1]} "+hours)
        end
      end

      # Minutes - assuming no one uses `1 minute`
      if ! match[2].nil? then
        duration_parts.push("#{match[2]} "+minutes)
      end

      # Hopefully no one uses seconds
      if ! match[3].nil? then
        duration_parts.push("#{match[3]} seconds")
      end

      return duration_parts.join(' ')
    end
  end
end

Liquid::Template.register_filter(Jekyll::DurationFilter)
