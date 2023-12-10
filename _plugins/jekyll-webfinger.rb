# frozen_string_literal: true

require './_plugins/gtn'

Jekyll::Hooks.register :site, :post_write do |site|
  # Make the directory
  Jekyll.logger.info 'Generating webfinger files'
  FileUtils.mkdir_p "#{site.dest}/api/fedi"

  Gtn::Contributors.list(site)
                   .select { |_k, v| v['fediverse'] }
                   .each do |k, v|
    # saving the outputs to
    # training-material/api/fedi/resource=acct%3Ahexylena%40galaxy.training.json

    subscribe_url = if v['fediverse_flavor'] == 'mastodon'
                      "#{v['fediverse'].gsub(%r{/@[a-z]*$}, '')}/authorize_interaction?uri={uri}"
                    else
                      "#{v['fediverse'].gsub(%r{/[a-z]*$}, '')}/ostatus_subscribe?acct={uri}"
                    end
    f2 = v['fediverse'].gsub(%r{https://(.*)/@?(.*)}, 'https://\1/users/\2')
    profile = {
      subject: "acct:#{k}@training.galaxyproject.org",
      aliases: [
        f2
      ],
      links: [
        {
          rel: 'http://webfinger.net/rel/profile-page',
          type: 'text/html',
          href: f2,
        },
        {
          rel: 'self',
          type: 'application/activity+json',
          href: f2 # ehhh
        },
        {
          rel: 'self',
          type: 'application/ld+json; profile="https://www.w3.org/ns/activitystreams"',
          href: f2,
        },
        {
          rel: 'http://ostatus.org/schema/1.0/subscribe',
          template: subscribe_url
        }
      ]
    }
    path = "#{site.dest}/api/fedi/resource=acct:#{k}@training.galaxyproject.org.json"
    File.write(path, profile.to_json)

    path = "#{site.dest}/api/fedi/resource=acct:@#{k}@training.galaxyproject.org.json"
    File.write(path, profile.to_json)
  end
end
