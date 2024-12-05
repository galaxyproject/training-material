# frozen_string_literal: true

OUR_PATH = __dir__
# two directories up
ROOT_PATH = File.expand_path(File.join(OUR_PATH, '..', '..'))

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
      Jekyll.logger.info '[GTN/Time/Mod] Filling Time Cache'
      cached_command
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

    def self.discover_caches
      # Really there should only be one, but maybe someone's been silly so
      # we'll just take the first one we find.
      Dir.glob('metadata/git-mod-*.txt').first
    end

    def self.generate_cache
      rev = `git rev-list -n 1 main`.strip

      if discover_caches.nil?
        File.write("metadata/git-mod-#{rev}.txt", command)
      else
        prev = discover_caches
        results = cached_command
        File.delete(prev)
        File.write("metadata/git-mod-#{rev}.txt", results)
      end
    end

    def self.cached_command
      return command if discover_caches.nil?

      Jekyll.logger.info '[GTN/Time/Mod] Using cached modification times'

      previous_commit = discover_caches.split('-').last.split('.').first
      previous = File.read(discover_caches)

      `git log --first-parent --name-only --pretty='GTN_GTN:%ct' #{previous_commit}..` + previous
    end

    def self.command
      `git log --first-parent --name-only --pretty='GTN_GTN:%ct'`
    end

    def self.time_cache
      @@TIME_CACHE
    end

    def self.commit_count_cache
      @@COMMIT_COUNT_CACHE
    end

    def self.clean_path(f)
      if f =~ %r{^\./}
        f[2..]
      elsif f =~ %r{^/}
        f.gsub(ROOT_PATH, '')
      else
        f
      end
    end

    def self.obtain_modification_count(f_unk)
      f = clean_path(f_unk)
      init_cache
      if @@COMMIT_COUNT_CACHE.key? f
        @@COMMIT_COUNT_CACHE[f]
      else
        0
      end
    end

    def self.obtain_time(f_unk)
      f = clean_path(f_unk)
      init_cache
      if @@TIME_CACHE.key? f
        @@TIME_CACHE[f]
      else
        begin
          # Non git file.
          @@TIME_CACHE[f] = File.mtime(f)
          Jekyll.logger.warn "[GTN/Time/Mod] No git cached time available for #{f}, defaulting to checkout"
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
    @@RENAMES = nil

    def self.chase_rename(path, depth: 0)
      if @@RENAMES.nil?
        self.init_cache
      end

      if @@RENAMES.key? path
        # TODO(hexylena)
        # This happens because it's the wrong datastructure, if there's a loop
        # in there, it'll just cycle through it endlessly.
        # This is obviously bad. But it'll do for now because it doesn't affect
        # any of our core files. We should replace this in the future.
        # This is why we have the grep_v below, to weed out the problematic files.
        if depth > 10
          Jekyll.logger.error "[GTN/Time/Pub] Too many renames for #{path}"
          path
        else
          chase_rename(@@RENAMES[path], depth: depth + 1)
        end
      else
        path
      end
    end

    def self.init_cache
      return unless @@TIME_CACHE.nil?

      @@TIME_CACHE = {}
      @@RENAMES = {}

      Jekyll.logger.info '[GTN/Time/Pub] Filling Publication Time Cache'
      cached_command
        .split('GTN_GTN:')
        .map { |x| x.split("\n\n") }
        .select { |x| x.length > 1 }
        .each do |date, files|
        files.split("\n").grep_v(/\.(png|json|_ga|jpg)/).each do |f|
          modification_type, path = f.split("\t")
          if modification_type == 'A'
            # Chase the renames.
            final_filename = chase_rename(path)
            @@TIME_CACHE[final_filename] = Time.at(date.to_i)
          elsif modification_type[0] == 'R'
            _, moved_from, moved_to = f.split("\t")
            @@RENAMES[moved_from] = moved_to # Point from the 'older' version to the newer.
          end
        end
      end
      # pp renames
    end

    def self.discover_caches
      # Really there should only be one, but maybe someone's been silly so
      # we'll just take the first one we find.
      Dir.glob('metadata/git-pub-*.txt').first
    end

    def self.generate_cache
      rev = `git rev-list -n 1 main`.strip

      if discover_caches.nil?
        File.write("metadata/git-pub-#{rev}.txt", command)
      else
        prev = discover_caches
        results = cached_command
        File.delete(prev)
        File.write("metadata/git-pub-#{rev}.txt", results)
      end
    end

    def self.cached_command
      return command if discover_caches.nil?

      Jekyll.logger.info '[GTN/Time/Pub] Using cached publication times'

      previous_commit = discover_caches.split('-').last.split('.').first
      previous = File.read(discover_caches)

      `git log --first-parent --name-status --diff-filter=AR --pretty='GTN_GTN:%ct' #{previous_commit}..` + previous
    end

    def self.command
      `git log --first-parent --name-status --diff-filter=AR --pretty='GTN_GTN:%ct' `
    end

    def self.time_cache
      @@TIME_CACHE
    end

    def self.clean_path(f)
      if f =~ %r{^\./}
        f[2..]
      elsif f =~ %r{^/}
        f.gsub(ROOT_PATH, '')
      else
        f
      end
    end

    def self.obtain_time(f_unk)
      f = clean_path(f_unk)
      init_cache
      if @@TIME_CACHE.key? f
        @@TIME_CACHE[f]
      else
        begin
          # Non git file.
          @@TIME_CACHE[f] = File.mtime(f)
          Jekyll.logger.warn "[GTN/Time/Pub] No git cached time available for #{f}, defaulting to checkout"
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

  puts ' Moved to bin/list-recently-modified.rb'
  # Gtn::PublicationTimes.init_cache
  # Gtn::PublicationTimes.time_cache.select do |_, v|
  #   # Things in last 6 months
  #   v > Time.now - (6 * 30 * 24 * 60 * 60)
  # end.map { |k, v| puts "#{v} #{k}" }
end
