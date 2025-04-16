require './_plugins/jekyll-topic-filter'
require 'jekyll'


module Gtn
  # Parse the git repo to get some facts
  module Hooks

    ##
    # Generate the by-tool pages
    # Params:
    # +site+:: Jekyll site object
    def self.by_tool(site)
      Jekyll.logger.debug "[GTN/Hooks/by_tool] Started"
      init_count = site.pages.size
      start_time = Time.now

      tools = Gtn::TopicFilter.list_materials_by_tool(site)
      tools.reject!{|tool, _| tool.include?('{{')}

      tools.each do |tool, tutorials|
        # tool: e.g. `saskia-hiltemann/krona_text/krona-text`

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

        # Redirect from the older, shorter IDs that have more potential for conflicts.
        if tool.include?('/')
          previous_id = tool.split('/')[0] + '/' + tool.split('/')[2]
        else
          previous_id = tool # No change
        end

        page2 = Jekyll::PageWithoutAFile.new(site, '', 'by-tool/', "#{tool.gsub('%20', ' ')}.html")
        page2.content = nil
        page2.data['layout'] = 'by_tool'
        page2.data['short_tool'] = tool
        page2.data['observed_tool_ids'] = ordered_tool_ids.map{|x| x[0..1]}.reverse
        page2.data['tutorial_list'] = tutorials['tutorials']
        page2.data['latest_tool_id'] = ordered_tool_ids.map{|x| x[0]}.last
        # page2.data['redirect_from'] = ["/by-tool/#{previous_id.gsub('%20', ' ')}"]
        site.pages << page2

        # TODO: For whatever reason the redirect_from does NOT work, even this
        # early in the hooks, so we're just going to write the file and call it
        # a day. Someone should fix this someday. My apologies for leaving it like this.
        if previous_id != tool
          page2 = Jekyll::PageWithoutAFile.new(site, '', 'by-tool/', "#{previous_id}.html")
          page2.content = nil
          page2.data['layout'] = 'by_tool'
          page2.data['short_tool'] = tool
          page2.data['observed_tool_ids'] = ordered_tool_ids.map{|x| x[0..1]}.reverse
          page2.data['tutorial_list'] = tutorials['tutorials']
          page2.data['latest_tool_id'] = ordered_tool_ids.map{|x| x[0]}.last
          site.pages << page2
        end

      end
      Jekyll.logger.info "[GTN/Hooks/by_tool] #{site.pages.size - init_count} pages added in #{Time.now - start_time}s"
    end
  end
end
