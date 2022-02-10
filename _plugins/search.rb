require 'json'
require 'liquid'
require './_plugins/colour-tags.rb'

module Jekyll
  class DumpSearchDataTag < Liquid::Tag

    def initialize(tag_name, text, tokens)
      super
      @text = text.strip
    end

    def getlist(tutorial, attr)
      if tutorial[attr].nil?
        []
      end

      tutorial.fetch(attr, [])
    end

    def render(context)
      site = context.registers[:site]
      topics = site.data.select{|k, v| v.is_a?(Hash) && v.has_key?('type')}
        .select{|k, v| ['use', 'admin-dev', 'basics'].include? v['type'] }

      results = {}
      topics.each{|k, topic|
        tutorials = site.data['cache_topic_filter'][k]
        tutorials.each{|tutorial|
          results[tutorial['url']] = {
            "topic" => topic['title'],
            "title" => tutorial['title'],
            "description" => tutorial['description'],
            "question" => getlist(tutorial, 'questions'),
            "objectives" => getlist(tutorial, 'objectives'),
            "tags" => getlist(tutorial, 'tags').map{|tag|
              %Q(<a class="label label-default" title="Show all tutorials tagged #{tag}" href="#{site.baseurl}/search?#{tag}" style="#{ColourTag.colour_tag tag}">#{tag}</a>)
            },
            "contributors" => getlist(tutorial, 'contributors').map{|c|
              site.data['contributors'].fetch(c, {}).fetch('name', c)
            }.join(', '),
            "level" => tutorial['level'],
            "time_estimation" => tutorial['time_estimation'],
            "url" => site.baseurl + tutorial['url'],
          }
        }
      }

      JSON.pretty_generate(results)
    end

  end
end

Liquid::Template.register_tag('dump_search_view', Jekyll::DumpSearchDataTag)
