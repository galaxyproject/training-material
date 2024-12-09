# frozen_string_literal: true

module Jekyll
  # This module contains a filter for converting a duration string into a human readable string.
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
      seconds = parse_rfc3339(duration)
      if seconds.nil?
        return duration
      end
      return fmt_duration(seconds)
    end

    def fmt_duration(seconds)
      d = resolve_hms(seconds)

      # Otherwise append english terms for the various parts
      duration_parts = []

      hour = 'hour'
      hours = 'hours'
      minutes = 'minutes'
      if @context.registers[:page]&.key?('lang') && (@context.registers[:page]['lang'] != 'en')
        lang = @context.registers[:page]['lang']
        hour = @context.registers[:site].data['lang'][lang]['hour']
        hours = @context.registers[:site].data['lang'][lang]['hours']
        minutes = @context.registers[:site].data['lang'][lang]['minutes']
      end

      # Hours
      if d[:hours] > 0
        if d[:hours] == 1
          duration_parts.push("#{d[:hours]} " + hour)
        else
          duration_parts.push("#{d[:hours]} " + hours)
        end
      end

      # Minutes - assuming no one uses `1 minute`
      duration_parts.push("#{d[:minutes]} " + minutes) if d[:minutes] > 0

      # Hopefully no one uses seconds
      duration_parts.push("#{d[:seconds]} seconds") if d[:seconds] > 0

      duration_parts.join(' ')
    end

    ##
    # Sum the durations correctly for multiple RFC3339 formatted durations.
    # Params:
    # +s+:: The RFC3339 formatted duration string
    # Returns:
    # +d+:: a number of seconds
    def parse_rfc3339(s)
      if s == 0
        return 0
      end

      # Match the different parts of the string, must match entire string or it
      # will fail.
      match = /^T?(?:([0-9]*)[Hh])*(?:([0-9]*)[Mm])*(?:([0-9.]*)[Ss])*$/.match(s)

      # If it doesn't match, pass through unedited so we don't cause unexpected
      # issues.
      if match.nil?
        Jekyll.logger.debug "[GTN/Durations]:", "Could not parse time: #{s}"
        return nil
      end

      return match[1].to_i * 3600 + match[2].to_i * 60 + match[3].to_i
    end

    ##
    # Turn a count of seconds into hours/minutes/seconds.
    # Params:
    # +Int+:: A number of seconds
    # Returns:
    # +Hash+:: A hash with keys for hours, minutes, and seconds
    #
    # Example:
    # resolve_hms(5400)
    # => { hours: 1, minutes: 30, seconds: 0 }
    def resolve_hms(seconds)
      # Normalize the total
      minutes = seconds / 60
      seconds = seconds % 60
      hours = minutes / 60
      minutes = minutes % 60

      { hours: hours, minutes: minutes, seconds: seconds }
    end

    ##
    # Sum the durations correctly for multiple RFC3339 formatted durations.
    # Params:
    # +materials+:: The GTN material objects
    # Returns:
    # +String+:: The human total duration
    def sum_duration(materials)
      Jekyll.logger.debug "[GTN/Durations]: sum durations with #{materials.length} materials."
      total = 0
      materials.each do |material|
        if ! material['time_estimation'].nil?
          Jekyll.logger.debug "    [GTN/Durations]: #{material['time_estimation']} #{material['title']} -> #{parse_rfc3339(material['time_estimation'])}"
          total += parse_rfc3339(material['time_estimation'])
        end
      end
      fmt_duration(total)
    end
  end
end

Liquid::Template.register_filter(Jekyll::DurationFilter)
