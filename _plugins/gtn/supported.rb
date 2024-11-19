# frozen_string_literal: true

require './_plugins/gtn/usegalaxy'

module Gtn
  # Handle tool support queries
  module Supported
    ##
    # Identify the servers that support a given tool list
    #
    # Params:
    # +data+:: The data from metadata/public-server-tools.json
    # +tool_list+:: The list of tools to check (either 'upload1' or
    #      'toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy10' style tools)
    # Returns:
    # +supported+:: A hash of supported servers, with the following structure:
    #  {
    #    exact: [server1, server2, ...],
    #    inexact: [server1, server2, ...]
    #  }
    def self.calculate(data, tool_list)
      # p "Calculating supported servers for this tool list"
      if data.nil? || data.empty? || tool_list.empty? || tool_list.nil?
        return {
          'exact' => [],
          'inexact' => Gtn::Usegalaxy.servers.map do |x|
                         x = x.transform_keys(&:to_s)
                         x['usegalaxy'] = true
                         x
                       end
        }
      end

      supported = { exact: {}, inexact: {} }
      tool_list.each do |tool|
        if tool.count('/') > 4
          # E.g. toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy10
          tool_id = tool.split('/')[0..4].join('/')
          tool_version = tool.split('/')[5..].join('/')
          # p "Checking #{tool_id} #{tool_version}... "

          if data['tools'].key?(tool_id)
            supported[:exact][tool] = data['tools'][tool_id][tool_version] if data['tools'][tool_id].key?(tool_version)

            supported[:inexact][tool] = data['tools'][tool_id].map { |_, v| v }.flatten.uniq.sort
          end
        elsif data['tools'].key?(tool)
          # E.g. 'upload1'
          # p "Checking #{tool}... "
          supported[:inexact][tool] = data['tools'][tool].map { |_, v| v }.flatten.uniq.sort
          supported[:exact][tool] = data['tools'][tool].map { |_, v| v }.flatten.uniq.sort
        end
      end

      # Exactly supporting servers:
      # this is the set of intersections across supported[:exact][*]
      exact_support = (0..data['servers'].length - 1).to_a
      tool_list.each do |tool|
        exact_support &= (supported[:exact][tool] || [])
      end

      # Inexactly supporting servers
      # Set of intersections across (union of supported[:exact] and supported[:inexact])
      inexact_support = (0..data['servers'].length - 1).to_a
      tool_list.each do |tool|
        et = supported[:exact][tool] || []
        it = supported[:inexact][tool] || []
        inexact_support &= (et | it)
      end
      # Remove the exactly supported ones because that would be extraneous we
      # check here if it's an array because the above code will occasionally
      # generate a 'false' value when merging sets.
      inexact_support -= exact_support

      usegalaxy_server_urls = Gtn::Usegalaxy.servers.map { |x| x[:url] }

      {
        'exact' => (exact_support || []).map do |id|
          data['servers'][id].update(
            { 'usegalaxy' => usegalaxy_server_urls.include?(data['servers'][id]['url']) }
          )
        end,
        'inexact' => (inexact_support || []).map do |id|
          data['servers'][id].update(
            { 'usegalaxy' => usegalaxy_server_urls.include?(data['servers'][id]['url']) }
          )
        end
      }
    end

    ##
    # Identify the servers that support a given tool list
    #
    # Params:
    # +data+:: The data from metadata/public-server-tools.json
    # +tool_list+:: The list of tools to check (either 'upload1' or
    #      'toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy10' style tools)
    # Returns:
    # +supported+:: A hash of supported servers, with the following structure:
    #  {
    #    servers: [a, b, c]
    #    tools: [
    #      {
    #        id: toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy10,
    #        version: 0.69.8+galaxy10,
    #        servers: [
    #          [0.69.8+galaxy10],
    #          [0.69.8+galaxy11, 0.69.8+galaxy12],
    #          nil
    #        ]
    #      }
    #    ]
    #  }
    def self.calculate_matrix(data, tool_list)
      structure = {
        'servers' => [],
        'tools' => [],
      }
      # p "Calculating supported servers for this tool list"
      if data.nil? || data.empty? || tool_list.empty? || tool_list.nil?
        return structure
      end

      structure['servers'] = data['servers']

      tool_list.each do |tool|
        tool_for_server = {
          'id' => tool,
          'servers' => []
        }

        if tool.count('/') > 4
          # E.g. toolshed.g2.bx.psu.edu/repos/iuc/circos/circos/0.69.8+galaxy10
          tool_id = tool.split('/')[0..4].join('/')
          tool_version = tool.split('/')[5..].join('/')
          tool_for_server['version'] = tool_version
        else
          tool_id = tool
          tool_version = 'local'
          tool_for_server['version'] = tool_version
        end

        data['servers'].each.with_index do |server, s_index|
          # If we've never seen this tool anywhere
          if ! data['tools'].key?(tool_id)
            res = {'state' => 'missing'}
          else
            # If we've seen this exact version on the current server/index
            if data['tools'][tool_id].key?(tool_version) && data['tools'][tool_id][tool_version].include?(s_index)
              res = {'state' => 'exact', 'version' => tool_version}
            elsif data['tools'][tool_id]
              res = {
                'state' => 'inexact', 
                'versions' => data['tools'][tool_id].select{|k, v| v.include?(s_index)}.keys
              }

              if res['versions'].length == 1 && res['versions'][0] == '_'
                res = {
                  'state' => 'local'
                }
              elsif res['versions'].length == 0
                res = {'state' => 'missing'}
              end
            end
          end

          res['server'] = server['url']
          tool_for_server['servers'].push(res)

        end

        structure['tools'].push(tool_for_server)
      end

      # TODO: would be nice to re-sort the servers by the number of exact matches
      structure
    end
  end
end

## Allow running from the CLI
if __FILE__ == $PROGRAM_NAME
  if ARGV.length.positive? && (ARGV[0] == 'test')
    require 'test/unit'
    # Testing for the class
    class IntersectionTest < Test::Unit::TestCase
      def test_exact
        data = {
          'servers' => { 0 => 's0', 1 => 's1', 2 => 's2' },
          'tools' => {
            'upload1' => {
              '_' => [0, 1, 2]
            },
            'ts/repos/hexy/lena/tool1' => {
              '1.0' => [0, 1],
              '2.0' => [0, 2]
            }
          }
        }

        # Inexact here is empty because all are exactly supporting
        assert_equal(Gtn::Supported.calculate(data, ['upload1']), { exact: %w[s0 s1 s2], inexact: [] })
        assert_equal(Gtn::Supported.calculate(data, ['upload1', 'ts/repos/hexy/lena/tool1/1.0']),
                     { exact: %w[s0 s1], inexact: %w[s2] })
        assert_equal(Gtn::Supported.calculate(data, ['upload1', 'ts/repos/hexy/lena/tool1/2.0']),
                     { exact: %w[s0 s2], inexact: %w[s1] })
        assert_equal(Gtn::Supported.calculate(data, %w[upload1 unknown-tool]),
                     { exact: %w[], inexact: %w[] })
      end
    end
  else
    server = ARGV[0]
    workflow = ARGV[1]

    def short_id(tool_id)
      if tool_id.count('/') > 4
        tool_id.split('/')[0..-2].join('/')
      else
        tool_id
      end
    end

    require 'json'
    tool_list = JSON.parse(`curl -s #{server}/api/tools?in_panel=False`).map{|x| x['id']}
    puts "Found #{tool_list.length} tools in #{server}"

    version_smash = tool_list.map{|x|
      if x.count('/') > 4
        [x.split('/')[0..-2].join('/'), x.split('/')[-1]]
      else
        [x, []]
      end
    }
    version_smash = version_smash.group_by{|x, y| x}.to_a.map{|k, v| [k, v.map{|x, y| y}]}.to_h

    workflow_tools = JSON.parse(File.open(workflow).read)['steps'].map{|_, x| x['tool_id']}
    workflow_tools.select!{|x| x && x.length.positive?}
    puts "Found #{workflow_tools.length} tools in the #{workflow}"

    workflow_tools.each do |tool|
      if tool_list.include?(tool)
        puts "✅ #{tool} is supported"
      else
        if version_smash.key?(short_id(tool))
          puts "❔ #{tool} is not supported, but #{version_smash[short_id(tool)]} are"
        else
          puts "❌ #{tool} is not supported"
        end
      end
    end

    metadata = {
      "servers" => [{"url" => server, "name" => "CLI Tested Server"}],
      "tools" => {}
    }
    version_smash.each do |k, v|
      if k.count('/') > 3
        metadata["tools"][k] = v.map{|x| [x, [0]]}.to_h
      else
        metadata["tools"][k] = {
          "_" => [0]
        }
      end
    end

    puts "GTN Rendering:"
    pp Gtn::Supported.calculate(metadata, workflow_tools)
  end
end
