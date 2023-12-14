# frozen_string_literal: true

module Gtn

  # Module for obtaining modification times of files.
  # It walks the git history to record the last time a file was modified.
  # This is faster than talking to the file system.
  module ModificationTimes
    @@TIME_CACHE = nil
    @@COMMIT_COUNT_CACHE = nil

    def self.init_cache
      return unless @@TIME_CACHE.nil?

      @@TIME_CACHE = {}
      @@COMMIT_COUNT_CACHE = Hash.new(0)
      puts '[GTN/Time/Mod] Filling Time Cache'
      `git log --name-only --pretty='GTN_GTN:%ct'`
        .split('GTN_GTN:')
        .map { |x| x.split("\n\n") }
        .select { |x| x.length > 1 }
        .each do |date, files|
        files.split(/\n/).each do |f|
          @@TIME_CACHE[f] = Time.at(date.to_i) if !@@TIME_CACHE.key? f
          @@COMMIT_COUNT_CACHE[f] += 1
        end
      end
    end

    def self.time_cache
      @@TIME_CACHE
    end

    def self.commit_count_cache
      @@COMMIT_COUNT_CACHE
    end

    def self.obtain_modification_count(f)
      init_cache
      if @@COMMIT_COUNT_CACHE.key? f
        @@COMMIT_COUNT_CACHE[f]
      else
        0
      end
    end

    def self.obtain_time(f)
      init_cache
      if @@TIME_CACHE.key? f
        @@TIME_CACHE[f]
      else
        begin
          # Non git file.
          @@TIME_CACHE[f] = File.mtime(f)
          @@TIME_CACHE[f]
        rescue StandardError
          Time.at(0)
        end
      end
    end
  end

  # Module for obtaining original publication times of files.
  # It walks the git history to record the last time a file was modified.
  # This is faster than talking to the file system.
  module PublicationTimes
    @@TIME_CACHE = nil

    def self.init_cache
      return unless @@TIME_CACHE.nil?

      @@TIME_CACHE = {}
      puts '[GTN/Time/Pub] Filling Publication Time Cache'
      `git log --first-parent --name-only --diff-filter=A --pretty='GTN_GTN:%ct' main`
        .split('GTN_GTN:')
        .map { |x| x.split("\n\n") }
        .select { |x| x.length > 1 }
        .each do |date, files|
        files.split(/\n/).select{|x| x =~ /\.(md|html)$/}.each do |f|
          @@TIME_CACHE[f] = Time.at(date.to_i) if !@@TIME_CACHE.key? f
        end
      end
    end

    def self.time_cache
      @@TIME_CACHE
    end

    def self.obtain_publication_time(f)
      init_cache
      if @@TIME_CACHE.key? f
        @@TIME_CACHE[f]
      else
        begin
          # Non git file.
          @@TIME_CACHE[f] = File.mtime(f)
          @@TIME_CACHE[f]
        rescue StandardError
          Time.at(0)
        end
      end
    end
  end

end

if $PROGRAM_NAME == __FILE__
  # Gtn::ModificationTimes.init_cache
  # pp Gtn::ModificationTimes.commit_count_cache

  Gtn::PublicationTimes.init_cache
  pp Gtn::PublicationTimes.time_cache
end
