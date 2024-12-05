# frozen_string_literal: true

require 'jekyll'
require 'time'

module Gtn
  # Parse the git repo to get some facts
  module Contributors
    @HAS_WARNED_ON = []

    ##
    # Returns contributors, regardless of whether they are 'contributor' or 'contributions' style
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def self.get_contributors(data)
      if data.key?('contributors')
        data['contributors']
      elsif data.key?('contributions')
        data['contributions'].values.flatten
      else
        []
      end
    end

    ##
    # Returns authors, only entities that are primary authors
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def self.get_authors(data)
      if data.key?('contributors')
        data['contributors'] || []
      elsif data.key?('contributions')
        data['contributions']['authorship'] || []
      else
        []
      end
    end

    ##
    # Returns event organisers
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def self.get_organisers(data)
      if data.key?('contributions') && data['contributions'].key?('organisers')
        data['contributions']['organisers']
      else
        []
      end
    end

    ##
    # Returns event instructors
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def self.get_instructors(data)
      if data.key?('contributions') && data['contributions'].key?('instructors')
        data['contributions']['instructors']
      else
        []
      end
    end

    ##
    # Get the non-author contributors of a material.
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def self.get_non_authors(material)
      if material.key?('contributors')
        []
      elsif material.key?('contributions')
        material['contributions']
          .reject { |k| k == 'funding' }
          .reject { |k| k == 'authorship' }
          .values.flatten.uniq
      else
        []
      end
    end

    ##
    # Get the funders (organisations) of a material.
    # Params:
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of contributor IDs
    def self.get_funders(site, data)
      self.get_all_funding(site, data).reject{ |f| site.data['grants'].key?(f) }
    end

    ##
    # Get the funders (grants) of a material.
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of grant IDs
    def self.get_grants(site, data)
      self.get_all_funding(site, data).select{ |f| site.data['grants'].key?(f) }
    end

    ##
    # Get the funders (grants + organisations) of a material.
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +data+:: +Hash+ of the YAML frontmatter from a material
    # Returns:
    # +Array+ of funder IDs mixing grants + organisations
    def self.get_all_funding(site, data)
      if data.key?('contributions') && data['contributions'].key?('funding')
        # The ones specifically in the Grants table
        data['contributions']['funding']
      else
        []
      end
    end

    # Convenience method to allow us to handle nil sites, and load directly
    # from disk ourselves.
    def self._load_file(site, category)
      if site.nil?
        Jekyll.logger.warn "[GTN/Contributor] Loading #{category} from disk, this access could be improved"
        File.open("_data/#{category}.yml", 'r') { |f| YAML.safe_load(f) }
      else
        site.data[category]
      end
    end

    ##
    # Map a contributor ID to their information and type
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +c+:: +String+ of contributor ID
    # Returns:
    # +Hash+ of contributor information
    # +String+ type of contributor (e.g. 'contributor', 'organisation', 'grant')
    def self.fetch(site, c, warn: false)
      if _load_file(site, 'contributors').key?(c)
        return ['contributor', site.data['contributors'][c]]
      elsif _load_file(site, 'organisations').key?(c)
        return ['organisation', site.data['organisations'][c]]
      elsif _load_file(site, 'grants').key?(c)
        return ['grant', site.data['grants'][c]]
      else
        if ! warn
          if ! @HAS_WARNED_ON.include?(c)
            Jekyll.logger.warn "Contributor #{c} not found"
            @HAS_WARNED_ON.push(c)
          end
        end
      end

      ['contributor', { 'name' => c }]
    end

    ##
    # Map a contributor ID to their information and type
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +c+:: +String+ of contributor ID
    # Returns:
    # +Hash+ of contributor information
    def self.fetch_contributor(site, c)
      fetch(site, c)[1]
    end

    ##
    # Map a contributor ID to their information and type
    # Params:
    # +site+:: +Jekyll::Site+ object
    # +c+:: +String+ of contributor ID
    # Returns:
    # +String+ of contributor name
    def self.fetch_name(site, c, warn: false)
      fetch(site, c, warn: warn)[1].fetch('name', c)
    end

    ##
    # List ALL contributors
    # Params:
    # +site+:: +Jekyll::Site+ object
    # Returns:
    # +Hash+ of contributors, grants, organisations merged together
    def self.list(site)
      site.data['contributors']
          .merge(site.data['grants'])
          .merge(site.data['organisations'])
          .reject { |c| c['halloffame'] == 'no' }
    end

    ##
    # Check if a specific contributor is a person or not
    # Params:
    # +c+:: +String+ of contributor ID
    # Returns:
    # +Boolean+ of whether the contributor is a contributor or not
    def self.person?(site, c)
      site.data['contributors'].key?(c)
    end

    ##
    # Check if a specific contributor is a grant or not
    # Params:
    # +c+:: +String+ of contributor ID
    # Returns:
    # +Boolean+ of whether the contributor is a grant or not
    def self.grant?(site, c)
      site.data['grants'].key?(c)
    end

    ##
    # Obtain the contributor's funding URL
    # Params:
    # +c+:: +String+ of contributor ID
    # Returns:
    # +Boolean+ of whether the contributor is a funder or not
    def self.fetch_funding_url(contributor)
      return contributor['funding_id'] if !contributor.key?('funding_database')

      case contributor['funding_database']
      when 'cordis'
        "https://cordis.europa.eu/project/id/#{contributor['funding_id']}"
      when 'erasmusplus'
        "https://erasmus-plus.ec.europa.eu/projects/search/details/#{contributor['funding_id']}"
      when 'ukri'
        "https://gtr.ukri.org/projects?ref=#{contributor['funding_id']}"
      when 'highergov'
        "https://www.highergov.com/contract/#{contributor['funding_id']}/"
      else
        Jekyll.logger.error "Unknown funding system #{contributor['funding_database']}. Please let us know so we can add support for it!"
        'ERROR'
      end
    end
  end
end
