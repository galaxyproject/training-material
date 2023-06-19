Jekyll::Hooks.register :site, :post_read do |site|
  site.config['javascript_bundles'].each do |name, resources|
    Jekyll.logger.info "Analysing JS Bundle #{name}"
    # Get the maximum last file modified time to use as the bundle timestamp
    bundle_timestamp = resources['resources'].map { |f| File.mtime(f).to_i }.max
    site.config['javascript_bundles'][name]['timestamp'] = bundle_timestamp
  end
end

Jekyll::Hooks.register :site, :post_write do |site|
  site.config['javascript_bundles'].each do |name, resources|
    bundle_path = "#{site.dest}/assets/js/bundle.#{name}.js"
    Jekyll.logger.info "Building JS bundle #{name} (#{bundle_path})"

    # Just concatenate them all together
    bundle = resources['resources'].map { |f| File.read(f) }.join("\n")
    # Write the bundle to the output directory
    File.write(bundle_path, bundle)
  end
end

module Jekyll
  # The main GTN function library
  module JsBundle
    def bundle_preloads(test)
      bundles = @context.registers[:site].config['javascript_bundles']
      baseurl = @context.registers[:site].config['baseurl']

      # Select the ones wishing to be preloaded
      bundles = bundles.select do |_name, bundle|
        bundle['preload'] == true
      end

      bundles.map do |name, bundle|
        bundle_path = "#{baseurl}/assets/js/bundle.#{name}.js?v=#{bundle['timestamp']}"

        "<link rel='preload' href='#{bundle_path}' as='script'>"
      end.join("\n")
    end

    def load_bundle(name)
      # TODO: does not load dev ever.
      if @context.registers[:site].config['environment'] == 'development'
        load_bundle_dev(name)
      else
        load_bundle_production(name)
      end
    end

    def load_bundle_dev(name)
      bundle = @context.registers[:site].config['javascript_bundles'][name]
      if bundle.nil?
        raise "Bundle #{name} not found in site config"
      end
      baseurl = @context.registers[:site].config['baseurl']

      bundle['resources'].map do |f|
        "<script src='#{baseurl}/#{f}'></script>"
      end.join("\n")
    end

    def load_bundle_production(name)
      bundle = @context.registers[:site].config['javascript_bundles'][name]
      if bundle.nil?
        raise "Bundle #{name} not found in site config"
      end
      baseurl = @context.registers[:site].config['baseurl']
      bundle_path = "#{baseurl}/assets/js/bundle.#{name}.js?v=#{bundle['timestamp']}"
      "<script src='#{bundle_path}'></script>"
    end
  end
end

Liquid::Template.register_filter(Jekyll::JsBundle)
