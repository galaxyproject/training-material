require 'json'


module Jekyll
  module JsonldFilter
    def to_jsonld(material, topic, site)
      gtn = {
        "@type": "Organization",
        "email": "#{site['email']}",
        "name": "Galaxy Training Network",
        "url": "https://galaxyproject.org/teach/gtn/"
      }
      if not topic then
        return '{}'
      end

      topic_desc = {
        "@type": "CreativeWork",
        "name": "#{topic['title']}",
        "description": "#{topic['summary']}",
        "url": "https://training.galaxyproject.org/#{site['baseurl']}/topics/#{topic['name']}/"
      }

      # aggregate everything
      data = {
        # Properties from Course
        "@context": "http://schema.org",
        "@type": "Course",
        #"courseCode" described below
        #"coursePrerequisites" described below
        #"educationalCredentialAwarded": ,
        #"hasCourseInstance": ,
        #"skillLevel" described below

        # Properties from CreativeWork
        #"about" described below
        "accessMode": ["textual", "visual"],
        "accessModeSufficient": ["textual", "visual"],
        #"accessibilityAPI": ,
        "accessibilityControl": ["fullKeyboardControl", "fullMouseControl"],
        "accessibilityFeature": ["alternativeText", "tableOfContents"],
        #"accessibilityHazard": [],
        "accessibilitySummary": "Short descriptions are present but long descriptions will be needed for non-visual users",
        #"accountablePerson":,
        #"aggregateRating":,
        #"alternativeHeadline":,
        #"associatedMedia":,
        "audience": {
          "@type": "EducationalAudience",
          "educationalRole": "students"
        },
        #"audio":,
        #"award":,
        #"author" described below
        #"character":,
        "citation": {
          "@type": "CreativeWork",
          "name": "Community-Driven Data Analysis Training for Biology",
          "url": "https://doi.org/10.1016/j.cels.2018.05.012"
        },
        #"comment":,
        #"commentCount":,
        #"contentLocation":,
        #"contentRating":,
        #"contentReferenceTime":,
        #"contributor" described below
        "copyrightHolder": gtn,
        #"copyrightYear":,
        #"correction":,
        #"creator":,
        #"dateCreated":,
        #"dateModified":,
        #"datePublished":,
        "discussionUrl": "#{site["gitter_url"]}",
        #"editor":,
        #"educationalAlignment":,
        #"educationalUse":,
        #"encoding":,
        #"encodingFormat":,
        #"exampleOfWork":,
        #"expires":,
        #"funder":,
        #"genre":,
        #"hasPart" described below
        "headline": "#{material['title']}",
        "inLanguage": {
            "@type": "Language",
            "name": "English",
            "alternateName": "en"
        },
        #"interactionStatistic":,
        "interactivityType": "mixed",
        "isAccessibleForFree": true,
        #"isBasedOn":,
        #"isFamilyFriendly":,
        #"isPartOf" described below
        #"keywords": described below
        #"learningResourceType" described below
        "license": "#{site['github_repository']}/blob/#{site['github_repository_branch']}/LICENSE.md",
        #"locationCreated":,
        #"mainEntity":,
        #"material":,
        #"mentions" described below
        #"offers":,
        #"position":,
        "producer": gtn,
        "provider": gtn,
        #"publication":,
        #"publisher":,
        #"publisherImprint":,
        #"publishingPrinciples":,
        #"recordedAt":,
        #"releasedEvent":,
        #"review":,
        #"schemaVersion":,
        #"sdDatePublished":,
        #"sdLicense":,
        #"sdPublisher":,
        "sourceOrganization": gtn,
        #"spatialCoverage":,
        #"sponsor":,
        #"temporalCoverage":,
        #"text":,
        #"thumbnailUrl":,
        #"timeRequired" described below
        #"translationOfWork":,
        #"translator": Google Translate???,
        #"typicalAgeRange":,
        #"version":,
        #"video":,
        #"workExample":,
        #"workTranslation":,

        # Properties from Thing
        #"additionalType":,
        #"alternateName":,
        #"description" described below
        #"disambiguatingDescription":,
        #"identifier":,
        #"image":,
        #"mainEntityOfPage":,
        #"name" described below
        #"potentialAction":,
        #"sameAs":,
        #"subjectOf":,
        # "url" described below
      }

      #info depending if tutorial, hands-on or slide level
      parts = []
      mentions = []
      description = []

      data['isPartOf'] = topic_desc

      if material['type'] == 'introduction' then
        data['courseCode'] = "#{material['topic_name']} / introduction / #{material['name']}"
        data['learningResourceType'] = "slides"
        data['name'] = "Introduction to '#{topic['title']}'"
        data['url'] = "https://training.galaxyproject.org/#{site['baseurl']}#{material['url']}"
        description.push("Slides for #{topic['title']}")
      elsif material['name'] == 'tutorial.md' or material['name'] == 'slides.html' then
        if material['name'] == 'tutorial.md' then
          data['courseCode'] = "#{material['topic_name']} / #{material['tutorial_name']} / hands-on"
          data['learningResourceType'] = "hands-on tutorial"
          data['name'] = "Hands-on for '#{material['title']}' tutorial"
        else
          data['courseCode'] = "#{material['topic_name']} / #{material['tutorial_name']} / slides"
          data['learningResourceType'] = "slides"
          data['name'] = "Slides for '#{material['title']}' tutorial"
        end
        data['url'] = "https://training.galaxyproject.org/#{site['baseurl']}#{material['url']}"

        # Time required
        if material.key?('time_estimation') and not material['time_estimation'].nil? then
          data['timeRequired'] = "PT#{material['time_estimation'].upcase}"
        end

        # Description with questions, objectives and keypoints
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

        # Keywords
        if material.key?('tags') then
          data['keywords'] = material['tags'].join(', ')
        end
        #Zenodo links
        if material.key?('zenodo_link') then
          mentions = mentions.push({
            "@type": "Thing",
            "url": "#{material['zenodo_link']}",
            "name": "Training data for #{material['title']} tutorial"
          })
        end
      end
      data['description'] = description.join('\n')

      # Course requirements (material + topic)
      reqs = []
      if topic.key?('requirements') then
        reqs.push(*topic['requirements'])
      end
      if material.key?('requirements') then
        reqs.push(*material['requirements'])
      end
      if !reqs.empty?
        coursePrerequisites = []
        for req in reqs do
          if req['type'] == "internal" then
            if req.key?('tutorials') then
              for tuto in req['tutorials'] do
                for page in site['pages'] do
                  if page['name'] == 'tutorial.md' or page['name'] == 'slides.html' then
                    if page['topic_name'] == req['topic_name'] and page['tutorial_name'] == tuto then
                      #slides
                      if page['name'] == 'slides.html' then
                        coursePrerequisites.push({
                          "@type": "Course",
                          "url": "https://training.galaxyproject.org/#{site['baseurl']}/topics/#{req['topic_name']}/tutorials/#{tuto}/slides.html",
                          "name": "#{page['title']}",
                          "description": "Slides for '#{page['title']}' tutorial",
                          "learningResourceType": "slides",
                          "interactivityType": "expositive",
                          "provider": gtn
                        })
                        if page['hands_on_url'] then
                          coursePrerequisites.push({
                            "@type": "Course",
                            "url": "#{page['hands_on_url']}",
                            "learningResourceType": "hands-on tutorial",
                            "interactivityType": "expositive",
                          })
                        end
                      end
                      #hands-on
                      if page['name'] == 'tutorial.md' then
                        coursePrerequisites.push({
                          "@type": "Course",
                          "url": "https://training.galaxyproject.org/#{site['baseurl']}/topics/#{req['topic_name']}/tutorials/#{tuto}/tutorial.html",
                          "name": "#{page['title']}",
                          "description": "Hands-on for '#{page['title']}' tutorial",
                          "learningResourceType": "hands-on tutorial",
                          "interactivityType": "expositive",
                          "provider": gtn
                        })
                      end
                    end
                  end
                end
              end
            else
              coursePrerequisites.push({
                "@type": "CreativeWork",
                "url": "https://training.galaxyproject.org/#{site['baseurl']}/topics/#{req['topic_name']}/",
                "name": "#{site['data'][req['topic_name']]['title']}",
                "description": "#{site['data'][req['topic_name']]['title']}",
                "provider": gtn
              })
            end
          elsif req['type'] == "external" then
            coursePrerequisites.push({
              "@type": "CreativeWork",
              "url": "#{req['link']}",
              "name": "#{req['title']}"
            })
          else
            coursePrerequisites.push("#{req['title']}")
          end
        end
        data['coursePrerequisites'] = coursePrerequisites
      end

      data['hasPart'] = parts

      # Add contributors/authors
      if material.key?('contributors') then
        contributors = material['contributors'].map{ |x|
          {
            "@type": "Person",
            "name": "#{site['data']['contributors'][x]['name']}" #expand to use real names
          }
        }
        data['author'] = contributors
        data['contributor'] = contributors
      end

      about = []
      about.push(topic_desc)
      if topic.key?('edam_ontology') then
        about.push({
          "@type": "DefinedTerm",
          "@id": "http://edamontology.org/#{topic['edam_ontology']}",
          "inDefinedTermSet": "http://edamontology.org",
          "termCode": "#{topic['edam_ontology']}",
          #"name": ,
          "url": "https://bioportal.bioontology.org/ontologies/EDAM/?p=classes&conceptid=http%3A%2F%2Fedamontology.org%2F#{topic['edam_ontology']}"
        })
      end
      data['about'] = about

      if material.key?('level') then
        data['skillLevel'] = material['level']
      end

      return JSON.pretty_generate(data)
    end
  end
end

Liquid::Template.register_filter(Jekyll::JsonldFilter)
