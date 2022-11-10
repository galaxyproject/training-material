require 'json'
require './_plugins/jekyll-topic-filter.rb'

module Jekyll
  class APIGenerator < Generator

    def get_contributors(data)
      if data.has_key?('contributors')
        return data['contributors']
      elsif data.has_key?('contributions')
        return data['contributions'].keys.map{|k| data['contributions'][k]}.flatten()
      end
    end

    def generate(site)

      # Full Bibliography
      puts "[GTN/API] Bibliography"
      page3 = PageWithoutAFile.new(site, "", "api/", "gtn.bib")
      page3.content = site.config['cached_global_bib'].to_s
      page3.data["layout"] = nil
      site.pages << page3

      def markdownify(site, text)
        site.find_converter_instance(
          Jekyll::Converters::Markdown
        ).convert(text.to_s)
      end

      def visitAndMarkdownify(site, f)
          if f.is_a?(Array)
              f.map!{|x| visitAndMarkdownify(site, x)}
          elsif f.is_a?(Hash)
              f = f.map{|k, v|
                  [k, visitAndMarkdownify(site, v)]
              }.to_h
          elsif f.is_a?(String)
              f = markdownify(site, f).strip().gsub(/<p>/, "").gsub(/<\/p>/, "")
          end
          f
      end

      def mapContributor(site, c)
        x = site.data['contributors']
          .fetch(c, {})
          .merge({"id" => c, "url" => site.config['url'] + site.config['baseurl'] + "/api/contributors/#{c}.json"})
        visitAndMarkdownify(site, x)
      end

      # Contributors
      puts "[GTN/API] Contributors"
      page2 = PageWithoutAFile.new(site, "", "api/", "contributors.json")
      page2.content = JSON.pretty_generate(site.data['contributors'].map{|c, _| mapContributor(site, c)})
      page2.data["layout"] = nil
      site.pages << page2
      site.data['contributors'].each{|c, _|
        page4 = PageWithoutAFile.new(site, "", "api/", "contributors/#{c}.json")
        page4.content = JSON.pretty_generate(mapContributor(site, c))
        page4.data["layout"] = nil
        site.pages << page4
      }

      # Trigger the topic cache to generate if it hasn't already
      puts "[GTN/API] Tutorials"
      TopicFilter.topic_filter(site, 'does-not-matter')
      TopicFilter.list_topics(site).map{|topic|

        out = site.data[topic].dup
        out['materials'] = TopicFilter.topic_filter(site, topic).map{|x|
          q = x.dup
          q['contributors'] = get_contributors(q).dup.map{|c| mapContributor(site, c)}

          q['urls'] = Hash.new

          if ! q['hands_on'].nil?
            q['urls']['hands_on'] = site.config['url'] + site.config['baseurl'] + "/api/topics/#{q['url'][7..-6]}.json"
          end

          if ! q['slides'].nil?
            q['urls']['slides'] = site.config['url'] + site.config['baseurl'] + "/api/topics/#{q['url'][7..-6]}.json"
          end

          # Write out the individual page
          page6 = PageWithoutAFile.new(site, "", "api/topics/", "#{q['url'][7..-6]}.json")
          page6.content = JSON.pretty_generate(q)
          page6.data["layout"] = nil
          site.pages << page6

          q
        }
        out['maintainers'] = out['maintainers'].map{|c| mapContributor(site, c)}

        page2 = PageWithoutAFile.new(site, "", "api/topics/", "#{topic}.json")
        page2.content = JSON.pretty_generate(out)
        page2.data["layout"] = nil
        site.pages << page2
      }

      topics = Hash.new
      puts "[GTN/API] Topics"
      # Individual Topic Indexes
      site.data.each_pair{|k, v|
        if v.is_a?(Hash) and v.has_key?('type') and v.has_key?('maintainers')

          topics[k] = {
            "name" => v['name'],
            "title" => v['title'],
            "summary" => v['summary'],
            "url" => site.config['url'] + site.config['baseurl'] + "/api/topics/#{k}.json",
            "maintainers" => v['maintainers'].map{|c| mapContributor(site, c) }
          }
        end
      }

      # Overall topic index
      page2 = PageWithoutAFile.new(site, "", "api/", "topics.json")
      page2.content = JSON.pretty_generate(topics)
      page2.data["layout"] = nil
      site.pages << page2

      def filterInteresting(layout)
        layout == 'tutorial_slides' or layout == 'base_slides' or layout == 'rdmbites_slides' or layout == 'tutorial_hands_on'
      end

      puts "[GTN/API] Tutorial and Slide pages"
      site.pages.select{|page| filterInteresting(page.data['layout']) }
        .each{|page|

        page5 = PageWithoutAFile.new(site, "", "api/topics/", "#{page.url[7..-6]}.json")
        p = page.data.dup
        p['contributors'] = get_contributors(p).dup.map{|c| mapContributor(site, c)}
        page5.content = JSON.pretty_generate(p)
        page5.data["layout"] = nil
        site.pages << page5
      }
      # Deploy the feedback file as well
      page2 = PageWithoutAFile.new(site, "", "api/", "feedback.json")
      page2.content = JSON.pretty_generate(site.data['feedback'])
      page2.data["layout"] = nil
      site.pages << page2

    end
  end
end

