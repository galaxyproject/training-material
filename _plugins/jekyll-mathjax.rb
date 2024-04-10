
Jekyll::Hooks.register :pages, :post_init do |page|
  page.data['js_requirements'] = {
    'mathjax' => page.content =~ /\$\$/,
    'mermaid' => page.content =~ /```mermaid/ || page.content =~ /pre class="mermaid"/ || page.data['layout'] == 'workflow-list',
  }

end
