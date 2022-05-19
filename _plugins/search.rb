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
      puts "[GTN/Search]"

      site = context.registers[:site]
      topics = site.data.select{|k, v| v.is_a?(Hash) && v.has_key?('type')}
        .select{|k, v| ['use', 'admin-dev', 'basics'].include? v['type'] }

      results = {}
      topics.each{|k, topic|
        tutorials = site.data['cache_topic_filter'][k]
        tutorials.each{|tutorial|
          results[tutorial['url']] = {
            "type" => 'Tutorial',
            "topic" => topic['title'],
            "title" => tutorial['title'],
            "contributors" => getlist(tutorial, 'contributors').map{|c|
              site.data['contributors'].fetch(c, {}).fetch('name', c)
            }.join(', '),
            "tags" => getlist(tutorial, 'tags').map{|tag|
              %Q(<a class="label label-default" title="Show all tutorials tagged #{tag}" href="#{site.baseurl}/search?query=#{tag}" style="#{ColourTag.colour_tag tag}">#{tag}</a>)
            },
            "url" => site.baseurl + tutorial['url'],
          }
        }
      }

      faqs = site.pages.select{|p| p.data['layout'] == 'faq'}
      faqs.each{|resource|
          results[resource['url']] = {
            "type" => 'FAQ',
            "topic" => 'FAQ',
            "title" => resource['title'],
            "contributors" => getlist(resource.data, 'contributors').map{|c|
              site.data['contributors'].fetch(c, {}).fetch('name', c)
            }.join(', '),
            "tags" => [],
            "url" => site.baseurl + resource['url'],
          }
      }

      JSON.pretty_generate(results)
    end

  end
end

Liquid::Template.register_tag('dump_search_view', Jekyll::DumpSearchDataTag)
