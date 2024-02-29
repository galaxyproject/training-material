require './_plugins/gtn'

# Write the commit log
Gtn::ModificationTimes.generate_cache
Gtn::PublicationTimes.generate_cache

# Gziped?
# require 'zlib'
# File.write("metadata/git-mod-#{rev}.txt.gz", Zlib.gzip(Gtn::ModificationTimes.command))
# File.write("metadata/git-pub-#{rev}.txt.gz", Zlib.gzip(Gtn::PublicationTimes.command))
# p Gtn::ModificationTimes.discover_caches
