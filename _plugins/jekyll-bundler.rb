require 'digest'

# While reading in all of the files for the site
# load the bundles and find out their last modified time
# so we can use that as a cache buster
Jekyll::Hooks.register :site, :post_read do |site|
  site.config['javascript_bundles'].each do |name, resources|
    # Get the maximum last file modified time to use as the bundle timestamp
    bundle_timestamp = resources['resources'].map { |f| File.mtime(f).to_i }.max
    site.config['javascript_bundles'][name]['timestamp'] = bundle_timestamp

    # This is inefficient since we read twice but it's also not that expensive.
    bundle = resources['resources'].map { |f| File.read(f) }.join("\n")
    hash = Digest::MD5.hexdigest(bundle)[0..7]
    site.config['javascript_bundles'][name]['hash'] = hash
    site.config['javascript_bundles'][name]['path'] = "/assets/js/bundle.#{hash}.js"

    Jekyll.logger.info "Analysing JS Bundle #{name} => #{bundle_timestamp} / #{hash}"
  end
end

# When writing the site, build the bundles
# It's basically "cat *.js > bundle.js"
# We don't need no fancy JS minification
# gzip probably does enough, everything else is pre-minified.
Jekyll::Hooks.register :site, :post_write do |site|
  site.config['javascript_bundles'].each do |name, resources|
    bundle_path = "#{site.dest}#{resources['path']}"
    Jekyll.logger.info "Building JS bundle #{name} => #{bundle_path}"

    # Just concatenate them all together
    bundle = resources['resources'].map { |f| File.read(f) }.join("\n")

    # Write the bundle to the output directory
    File.write(bundle_path, bundle)
  end
end

module Jekyll
  # The main GTN function library
  module JsBundle
    # Return the preloads for the bundles, when in production
    # +test+:: ignore this
    # Returns the HTML to load the bundle
    #
    # Example:
    # {{ 'load' | bundle_preloads }}
    def bundle_preloads(_test)
      if Jekyll.env == 'production'
        bundle_preloads_prod
      else
        ''
      end
    end

    # (Internal) Return the production preloads for the bundles
    def bundle_preloads_prod
      bundles = @context.registers[:site].config['javascript_bundles']
      baseurl = @context.registers[:site].config['baseurl']

      # Select the ones wishing to be preloaded
      bundles = bundles.select do |_name, bundle|
        bundle['preload'] == true
      end

      bundles.map do |_name, bundle|
        bundle_path = "#{baseurl}#{bundle['path']}?v=#{bundle['timestamp']}"
        "<link rel='preload' href='#{bundle_path}' as='script'>"
      end.join("\n")
    end

    # Load a specific bundle, in liquid
    # +name+:: the name of the bundle to load
    # Returns the HTML to load the bundle
    #
    # Example:
    # {{ 'main' | load_bundle }}
    def load_bundle(name)
      if Jekyll.env == 'development'
        load_bundle_dev(name)
      else
        load_bundle_production(name)
      end
    end

    def load_bundle_dev(name)
      bundle = @context.registers[:site].config['javascript_bundles'][name]
      raise "Bundle #{name} not found in site config" if bundle.nil?

      baseurl = @context.registers[:site].config['baseurl']

      bundle['resources'].map do |f|
        "<script src='#{baseurl}/#{f}'></script>"
      end.join("\n")
    end

    def load_bundle_production(name)
      bundle = @context.registers[:site].config['javascript_bundles'][name]
      raise "Bundle #{name} not found in site config" if bundle.nil?

      baseurl = @context.registers[:site].config['baseurl']
      bundle_path = "#{baseurl}#{bundle['path']}?v=#{bundle['timestamp']}"
      "<script src='#{bundle_path}'></script>"
    end
  end
end

Liquid::Template.register_filter(Jekyll::JsBundle)
