module Gtn
  module ModificationTimes
    @@TIME_CACHE = nil

    def self.init_cache
      return unless @@TIME_CACHE.nil?

      @@TIME_CACHE = {}
      puts '[GTN/MOD] Filling Time Cache'
      results = `git log --name-only --pretty='GTN_GTN:%ct'`.split('GTN_GTN:')
      results.map! { |x| x.split(/\n\n/) }
      results.select! { |x| x.length > 1 }
      results.each do |date, files|
        files.split(/\n/).each do |f|
          @@TIME_CACHE[f] = Time.at(date.to_i) if !@@TIME_CACHE.has_key? f
        end
      end
    end

    def self.time_cache
      @@TIME_CACHE
    end

    def self.obtain_time(f)
      init_cache
      if @@TIME_CACHE.has_key? f
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

if $0 == __FILE__
  Gtn::ModificationTimes.init_cache
  pp Gtn::ModificationTimes.time_cache
end
