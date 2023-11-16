def jupyter_pre_render(site)
  nil
end

def jupyter_post_write(site)
  nil
end

module Jekyll
  class RmarkdownGenerator < Generator
    def generate(site)
      puts "Notebooks disabled"
    end
  end

  class JupyterNotebookGenerator < Generator
    def generate(site)
      puts "Notebooks disabled"
    end
  end
end
