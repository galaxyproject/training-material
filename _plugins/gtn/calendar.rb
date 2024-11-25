require 'icalendar'
require 'kramdown'
require './_plugins/util'

module Gtn
  # Parse calendars into GTN style event objects
  module Calendar

    def self.parse(calendar, after: nil)
      # Open a file or pass a string to the parser
      contents = request(calendar).body

      # Parser returns an array of calendars because a single file
      # can have multiple calendars.
      cals = Icalendar::Calendar.parse(contents)
      cal = cals.first

      events = cal.events
      puts "after: #{after}"
      if after
        cal.events.select! { |e| e.dtstart.strftime("%Y-%m-%d") > after.strftime("%Y-%m-%d") }
      end

      events.map do |e|
        self.extract_event(e)
      end
    end

    def self.extract_event(e)
      e.description.force_encoding('UTF-8')
      e.summary.force_encoding('UTF-8')
      {
        'id' => e.summary.downcase.gsub(/\s+/, '-').gsub(/[^a-z0-9-]/, ''),
        'date' => e.dtstart,
        'meta' => {
          'layout' => 'event',
          'title' => e.summary.to_s,
          'date_start' => e.dtstart.to_date,
          'date_end' => e.dtend.to_date,
          'tags' => ['google-calendar'],
          'mode' => 'online',
          'async' => false,
        },
        'content' => Kramdown::Document.new(e.description, :html_to_native => true)
          .to_kramdown
          .gsub('</p>', "\n\n")
          .gsub('<p>', "\n\n")
          .gsub('/^\s*$/','')
          .gsub(/\n{3,}/, "\n\n")
          .strip
      }
    end
  end
end
