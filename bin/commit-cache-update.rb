require 'jekyll'
require './_plugins/gtn/mod'

# Write the commit log
Gtn::ModificationTimes.generate_cache
Gtn::PublicationTimes.generate_cache
