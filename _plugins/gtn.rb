require './_plugins/gtn/boxify'
require './_plugins/gtn/mod'
require './_plugins/gtn/metrics'


module Jekyll
  module GtnFunctions

    # These two could be unified tbh
    def last_modified_at(page)
      Gtn::ModificationTimes.obtain_time(page['path'])
    end

    def gtn_mod_date(path)
      # Automatically strips any leading slashes.
      Gtn::ModificationTimes.obtain_time(path.gsub(/^\//, ''))
    end

    def how_many_topic_feedbacks(feedback, name)
      feedback.select{|x| x["topic"] == name}.length
    end

    def how_many_tutorial_feedbacks(feedback, name)
      feedback.select{|x| x["tutorial"] == name}.length
    end

    def filter_authors(contributors, contributions)
      if not contributors.nil?
        return contributors
      else
        return contributions["authorship"]
      end
    end

    def get_default_link(material)
      url = nil

      if material['type'] == "introduction"
        subfolder = 'slides'
      else
        subfolder = 'tutorials'
      end

      if material['slides']
        url = "topics/#{material['topic_name']}/#{subfolder}/#{material['tutorial_name']}"
        if material['type'] != "introduction"
          url += "/slides.html"
        else
          url += ".html"
        end
      end

      if material['hands_on']
        if material['hands_on'] != "external" && material['hands_on'] != ""
          url = "topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/tutorial.html"
        end
      end

      url
    end

  end
end

Liquid::Template.register_filter(Jekyll::GtnFunctions)

if $0 == __FILE__
  result = Gtn::ModificationTimes.obtain_time(ARGV[0].gsub(/^\//, ''))
  puts "Modification time of #{ARGV[0].gsub(/^\//, '')} is #{result}"
end
