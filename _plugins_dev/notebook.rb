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
