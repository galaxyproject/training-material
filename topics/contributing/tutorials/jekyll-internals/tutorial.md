---
layout: tutorial_hands_on

title: "Jekyll Internals (Notes)"
questions:
  - "What internals of Jekyll should I know about"
objectives:
  - "Retain hope in the face of Ruby"
time_estimation: "2h"
key_points:
  - Jekyll is complicated but lets you do SO many nice things easily.
  - Please make more plugins.
contributors:
  - hexylena
---

This will be @hexylena's assorted notes on how Jekyll works internally to help other GTN maintainers.

## Hooks

[Hooks](https://jekyllrb.com/docs/plugins/hooks/) are super useful, we use them to decrease build time in a couple places. If you want to write content 'by hand' that escapes Liquid's processing, then you'll be interested in Hooks.

E.g. we write out some API endpoints, these absolutely do not need to be processed by liquid and the kramdown pipeline. It'd be a waste of time. So you can write them directly, after the rest of the site has been written:

```ruby
Jekyll::Hooks.register :site, :post_write do |site|
  # Make the directory
  Jekyll.logger.info 'Generating API files'
  FileUtils.mkdir_p "#{site.dest}/api/example"
  profile = {hello: "world"}

  path = "#{site.dest}/api/example/test.json"
  File.write(path, profile.to_json)
  end
end
```

> <tip-title>Inclusion in sitemap.xml</tip-title>
> Note that these files will not be included in the sitemap, which is potentially sub-optimal! That's generated earlier and probably looks at `pages` and `posts`, so you'd require that it's registered with one of those two. Perhaps you could move this up slightly in the pipeline to avoid that?
{: .tip}

too easy!

That's a `post_write` hook, but there's a wealth of other ones. E.g. if you want to munge a page before it's written, so you don't have to figure out where it should be written to, you can do a `:post_render` hook of `:pages`

```ruby
Jekyll::Hooks.register :pages, :post_render do |page|
  if page.output =~ /APPLE/
    page.output = page.output.gsub('APPLE', 'BANANA')
  end
end
```

This would apply to every single page that has the word APPLE in all caps in the page's text, but you could use this for other things! We use this to fix any remaining boxes before the html is written to disk.

There's also `:pre_render` hooks if you want to muck about with the markdown before it's templated out:

```
Jekyll::Hooks.register :posts, :pre_render do |post, _out|
  post.data['author'] = get_authors(post.data).map { |c| lookup_name(c, post.site) }.join(', ')
  post.data['image'] = post.data['cover']
end
```

Here we map all `post` objects (not pages, specifically the RSS posts) to copy over a flattened version of the author's actual names from the internal complicated structure we use for authorship. Same for the post's cover, it's duplicated into `image`.

This is because our RSS plugin (`jekyll-feed`) expects those names and provides little configurability, so we just do a quick lil' addition of new fields directly onto the post object.

## Plugins

They're just a ruby file, you can do whatever you want in them. Generally you'll want either a Jekyll hook, or some bit of code that should be run. The code will be run essentially at a random point in the array of plugins, so if you need something pre-processed, modify that plugin to ensure ordering, or, yeah.

## Please Do This

Do | Do not (be like Helena)
--- | ---
`Jekyll.logger.debug` (warn, info, error) | `puts `
