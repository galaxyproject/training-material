Jekyll::Hooks.register :site, :post_write do |site|
  # Make the directory
  Jekyll.logger.info "Generating webfinger files"
  FileUtils.mkdir_p "#{site.dest}/api/fedi"

  site.data['contributors']
    .select { |k, v| v['fediverse'] }
    .each do |k, v|
      # saving the outputs to
      # training-material/api/fedi/resource=acct%3Ahexylena%40galaxy.training.json
      path = "#{site.dest}/api/fedi/resource=acct:#{k}@training.galaxyproject.org.json"

      subscribe_url = if v['fediverse_flavor'] == 'mastodon'
                        "https://#{v['fediverse'].split('@')[1]}/users/#{k}/outbox"
                        v['fediverse'].gsub(/\/@[a-z]*$/, '') + "/authorize_interaction?uri={uri}"
                      else
                        v['fediverse'].gsub(/\/[a-z]*$/, '') + "/ostatus_subscribe?acct={uri}"
                      end
      profile = {
        "subject": "acct:#{k}@training.galaxyproject.org",
        "aliases": [
          v['fediverse'],
        ],
        "links": [
          {
            "rel": "http://webfinger.net/rel/profile-page",
            "type": "text/html",
            "href": v['fediverse']
          },
          {
            "rel": "self",
            "type": "application/activity+json",
            "href": v['fediverse'] #ehhh
          },
          {
            "rel": "self",
            "type": "application/ld+json; profile=\"https://www.w3.org/ns/activitystreams\"",
            "href": v['fediverse']
          },
          {
            "rel": "http://ostatus.org/schema/1.0/subscribe",
            "template": subscribe_url
          }
        ]
      }
      File.write(path, profile.to_json)
  end
end
