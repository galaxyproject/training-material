require './_plugins/jekyll-topic-filter'
require 'jekyll'


module Gtn
  # Parse the git repo to get some facts
  module Hooks

    def self.by_tool(site)
      Jekyll.logger.info "[GTN/Hooks/by_tool] Started"
      TopicFilter.list_materials_by_tool(site).each do |tool, tutorials|
        page2 = Jekyll::PageWithoutAFile.new(site, '', 'by-tool/', "#{tool.gsub('%20', ' ')}.html")
        page2.content = nil
        page2.data['layout'] = 'by_tool'
        page2.data['short_tool'] = tool

        ordered_tool_ids = tutorials['tool_id']
          .map{|x| 
            if x[0] == x[1]
              # TODO: collect versions of builtins.
              [x[0], '0.0.0'] # Fake version for local only tools
            else
              x
            end
          }
          .reject{|x| x[0] == x[1]}
          .map{|x| [x[0], x[1], Gem::Version.new(fix_version(x[1]))]}
          .sort_by{|x| x[2]}

        page2.data['observed_tool_ids'] = ordered_tool_ids.map{|x| x[0..1]}.reverse
        page2.data['tutorial_list'] = tutorials['tutorials']
        page2.data['latest_tool_id'] = ordered_tool_ids.map{|x| x[0]}.last

        # Redirect from the older, shorter IDs that have more potential for conflicts.
        if tool.include?('/')
          previous_id = tool.split('/')[0] + '/' + tool.split('/')[2]
        else
          previous_id = tool # No change
        end
        page2.data['redirect_from'] = ["/by-tool/#{previous_id.gsub('%20', ' ')}"]
        site.pages << page2
      end

    end
  end
end
