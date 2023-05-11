require 'jekyll'
require 'time'

module Gtn
  module Git
    def self.cache
      @@cache ||= Jekyll::Cache.new('Git')
    end

    ##
    # Discover git-related facts and ensure they're cached
    # Params:
    # +site+:: The +Jekyll::Site+ object
    # Returns:
    # +Hash+:: A hash of git facts like the current revision, tags, etc.
    def self.discover
      self.cache.getset('facts') do
        self._discover
      end
    end

    def self._discover
      # Cache the git facts

      begin
        git_head = File.open(File.join('.git', 'HEAD')).read.strip.split(' ')[1]
        git_head_ref = File.open(File.join('.git', git_head)).read.strip
      rescue
        git_head_ref = 'none'
      end

      begin
        tags = `git tag -l`.strip.split.sort
      rescue
        tags = []
      end

      first_commit = Date.parse('2015-06-29')
      today = Date.today

      {
        'git_revision' => git_head_ref,
        'git_revision_short' => git_head_ref[0..6],
        'gtn_fork' => ENV.fetch('GTN_FORK', 'galaxyproject'),
        'git_tags' => tags,
        'git_tags_recent' => tags.reverse[0..2],
        'age' => (today - first_commit).to_f / 365.25,
      }
    end
  end
end
