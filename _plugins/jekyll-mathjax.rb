Jekyll::Hooks.register :pages, :post_init do |page|
  page.data['mathjax'] = page.content =~ /\$\$/
end
