require 'json'


module Jekyll
  module JsonldFilter
    def to_jsonld(material)
      data = {
        "@context": "http://schema.org",
        "@type": "CreativeWork",
        "name": material['title'],
        "accessModeSufficient": "visual",
        "accessibilityControl": "fullMouseControl",
        "accessibilityFeature": ["alternativeText", "tableOfContents"],
        "audience": {
            "@type": "EducationalAudience",
            "educationalRole": "students"
        },
        "inLanguage": {
            "@type": "Language",
            "name": "English",
            "alternateName": "en"
        },
        "interactivityType": "mixed",
        "isAccessibleForFree": true,
        # TODO
        "license": "https://github.com/galaxyproject/training-material/blob/master/LICENSE.md",
        "learningResourceType": "#{material['type']}",
        "url": "#{material['url']}",
      }

      # Zenodo links
      mentions = []
      if material.key?('zenodo_link') then
        mentions = mentions.push({
          "@type": "Thing",
          "url": "#{material['zenodo_link']}",
          "name": "Training data for #{material['title']} tutorial"
        })
      end

      # Workflows
      if material.key?('workflows') then
        mentions = mentions.push({
          "@type": "Thing",
          # TODO
          "url": "https://github.com/galaxyproject/training-material/tree/master/topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/workflows/",
          "name": "Workflow for #{material['title']} tutorial"
        })
      end
      data['mentions'] = mentions

      if material.key?('time_estimation') and not material['time_estimation'].nil? then
        data['timeRequired'] = "PT#{material['time_estimation'].upcase}"
      end

      # Add contributors
      if material.key?('contributors') then
        data['contributor'] = material['contributors'].map{ |x|
          {
            "@type": "Person",
            "name": "#{x}"
          }
        }
      end

      # Keywords
      if material.key?('tags') then
        data['keywords'] = material['tags'].join(', ')
      end

      # Description
      description = []
      if material.key?('questions') and not material['questions'].nil? and material['questions'].length > 0 then
        questions = material['questions'].join("\n - ")
        description.push("The questions this #{material['type']} addresses are:\n - #{questions}\n\n")
      end

      if material.key?('objectives') and not material['objectives'].nil? and material['objectives'].length > 0 then
        objectives = material['objectives'].join("\n - ")
        description.push("The objectives are:\n - #{objectives}\n\n")
      end

      if material.key?('keypoints') and not material['keypoints'].nil? and material['keypoints'].length > 0 then
        keypoints = material['keypoints'].join("\n - ")
        description.push("The keypoints are:\n - #{keypoints}\n\n")
      end
      if description.length > 0 then
        data['description'] = description.join('\n')
      end

      # various subtypes
      if material['type'] == 'tutorial' then
        if material.key?('hands_on') then
          parts = []
          if material.key?('slides') then
            slide_part = {
              "@type": "CreativeWork",
              "url": "https://github.com/galaxyproject/training-material/tree/master/topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/tours/",
              "name": "Slides for #{material['title']}",
              "learningResourceType": "slides",
              "interactivityType": "expositive"
            }
            parts.push(slide_part)
          end

          if material.key?('galaxy_tour') then
            tour_part = {
              "@type": "CreativeWork",
              "url": "https://github.com/galaxyproject/training-material/tree/master/topics/#{material['topic_name']}/tutorials/#{material['tutorial_name']}/tours/",
              "name": "Galaxy Tour for #{material['title']}",
              "interactivityType": "active",
              "learningResourceType": "interactive-tour"
            }
            parts.push(tour_part)
          end

          data['hasPart'] = parts
        end
      end


      return JSON.pretty_generate(data)
    end
  end
end

Liquid::Template.register_filter(Jekyll::JsonldFilter)
