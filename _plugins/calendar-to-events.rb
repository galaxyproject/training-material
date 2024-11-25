# frozen_string_literal: true

require 'yaml'

require './_plugins/gtn'

# Basically like `PageWithoutAFile`, we just write out the ones we'd created earlier.
Jekyll::Hooks.register :site, :post_write do |site|
  site.config['calendars'].each do |_, calendar|
    url = calendar['url']
    tags = calendar['tags']
    organisers = calendar['organisers']
    Gtn::Calendar.parse(url, after: calendar['after']).each do |event|
      # Write out the event as a JSON file
      path = "events/#{event['date'].strftime('%Y-%m-%d')}-#{event['id']}.md"

      File.open(path, 'w') do |file|
        event['meta']['tags'].push(tags).flatten!
        event['meta']['contributions'] = {'organisers' => organisers}

        file.write(event['meta'].to_yaml)
        file.write("---\n")
        file.write(event['content'])
      end
    end
  end
end
