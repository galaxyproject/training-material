module Jekyll
  module DurationFilter

    ##
    # This function converts a duration string into a human readable string.
    # Params:
    # +duration+:: The duration string to convert (e.g. 1H30M, RFC 3339 formatted minus the leading P/T)
    # Returns:
    # +String+:: The human readable string
    #
    # Example:
    #  {{ "T1H30M" | duration_to_human }}
    #  => "1 hour 30 minutes"
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
