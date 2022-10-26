module Gtn
  module ModificationTimes
    @@TIME_CACHE = nil

    def self.init_cache
      if @@TIME_CACHE.nil?
        @@TIME_CACHE = Hash.new
        puts "[GTN/MOD] Filling Time Cache"
        results = `git log --name-only --pretty='GTN_GTN:%ct'`.split('GTN_GTN:')
        results.map!{|x| x.split(/\n\n/)}
        results.select!{|x| x.length > 1}
        results.each{|date, files|
          files.split(/\n/).each {|f|
            if not @@TIME_CACHE.has_key? f
              @@TIME_CACHE[f] = Time.at(date.to_i)
            end
          }
        }
      end
    end

    def self.obtain_time(f)
      self.init_cache
      if @@TIME_CACHE.has_key? f
        @@TIME_CACHE[f]
      else
        begin
          # Non git file.
          @@TIME_CACHE[f] = File.mtime(f)
          @@TIME_CACHE[f]
        rescue
          Time.at(0)
        end
      end
    end
  end
end
