require './_plugins/gtn/boxify'
require './_plugins/gtn/mod'


module Jekyll
  module GtnFunctions
    def gtn_mod_date(path)
      Gtn::ModificationTimes.obtain_time(path)
    end
  end
end

Liquid::Template.register_filter(Jekyll::GtnFunctions)
