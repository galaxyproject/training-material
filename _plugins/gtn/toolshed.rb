# frozen_string_literal: true

module Gtn
  # Handle toolshed yaml formatting for ephemeris
  module Toolshed
    ##
    # Prepare a tool list for installation
    #
    # Params:
    # +data+:: The data from metadata/toolshed-revisions.json
    # +tool_list+:: The list of tools to check (either 'upload1' or
    #               'toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy10' style tools)
    # +topic+:: The topic to install the tools under
    # Returns:
    # +supported+:: A string of the admin install, ready for ephemeris
    def self.format_admin_install(data, tool_list, topic, tool_cats)
      # p "Calculating supported servers for this tool list"
      return {} if data.nil? || data.empty?

      tools = tool_list.select { |t| data.key? t }.map do |tool|
        tool_info = data[tool]
        {
          'name' => tool_info[1],
          'owner' => tool_info[0],
          'revisions' => tool_info[2],
          'tool_panel_section_label' => tool_cats["#{tool_info[0]}/#{tool_info[1]}"] || topic,
          'tool_shed_url' => 'https://toolshed.g2.bx.psu.edu/',
        }
      end

      {
        'install_tool_dependencies' => true,
        'install_repository_dependencies' => true,
        'install_resolver_dependencies' => true,
        'tools' => tools
      }
    end
  end
end
