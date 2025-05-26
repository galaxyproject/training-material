#!/usr/bin/env ruby
require './_plugins/gtn'

months = (ARGV[1] || 6).to_i

Gtn::PublicationTimes.init_cache
Gtn::PublicationTimes.time_cache.select do |_, v|
  # Things in last 6 months
  v > Time.now - (months * 30 * 24 * 60 * 60)
end.map { |k, v| puts "#{v} #{k}" }
