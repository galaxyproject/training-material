# frozen_string_literal: true

module Jekyll
  module Tags
    class SnippetIncludeTag < IncludeTag

      def markdownify(text)
        @site.find_converter_instance(
          Jekyll::Converters::Markdown
        ).convert(text.to_s)
      end

      def render(context)

        @site ||= context.registers[:site]

        file = render_variable(context) || @file
        validate_file_name(file)

        @site.inclusions[file] ||= locate_include_file(file)
        inclusion = @site.inclusions[file]

        add_include_to_dependency(inclusion, context) if @site.config["incremental"]

        context.stack do
          context["include"] = parse_params(context) if @params
          x = inclusion.render(context)
          x.gsub!(/\A---(.|\n)*?---/, '')
          markdownify(x).gsub(/\r/, '').gsub(/\n/, '')
        end
      end

      private

      def locate_include_file(file)
        @site.includes_load_paths.each do |dir|
          path = PathManager.join(dir, file)
          return Inclusion.new(@site, dir, file) if valid_include_file?(path, dir)
        end
        raise IOError, could_not_locate_message(file, @site.includes_load_paths, @site.safe)
      end

      def valid_include_file?(path, dir)
        File.file?(path) && !outside_scope?(path, dir)
      end

      def outside_scope?(path, dir)
        @site.safe && !realpath_prefixed_with?(path, dir)
      end

      def realpath_prefixed_with?(path, dir)
        File.realpath(path).start_with?(dir)
      rescue StandardError
        false
      end

      def add_include_to_dependency(inclusion, context)
        return unless context.registers[:page]&.key?("path")

        @site.regenerator.add_dependency(
          @site.in_source_dir(context.registers[:page]["path"]),
          inclusion.path
        )
      end
    end
  end
end

Liquid::Template.register_tag("snippet", Jekyll::Tags::SnippetIncludeTag)
