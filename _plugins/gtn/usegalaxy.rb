module Gtn
  # Data about the UseGalaxy.* servers
  module Usegalaxy
    @@SERVERS = [
      { name: 'UseGalaxy.eu', url: 'https://usegalaxy.eu', id: 'eu', human: 'Galaxy Europe' },
      { name: 'UseGalaxy.org', url: 'https://usegalaxy.org', id: 'us', human: 'Galaxy Main' },
      { name: 'UseGalaxy.org.au', url: 'https://usegalaxy.org.au', id: 'au', human: 'Galaxy Australia' },
      { name: 'UseGalaxy.fr', url: 'https://usegalaxy.fr', id: 'fr', human: 'Galaxy France' }
    ]

    def self.servers
      @@SERVERS
    end
  end
end
