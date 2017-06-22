## ========================================
## Commands to build websites and slides

# Settings
MAKEFILES=Makefile $(wildcard *.mk)
JEKYLL=bundle exec jekyll
PARSER=bin/markdown_ast.rb
DST=_site

# Controls
.PHONY : commands clean files
all : commands

## commands         : show all commands.
commands :
	@grep -h -E '^##' ${MAKEFILES} | sed -e 's/## //g'

## serve            : run a local server.
serve :
	${JEKYLL} serve

## site             : build files but do not run a server.
site :
	${JEKYLL} build

## clean            : clean up junk files.
clean :
	@rm -rf ${DST}
	@rm -rf .sass-cache
	@find . -name .DS_Store -exec rm {} \;
	@find . -name '*~' -exec rm {} \;

## install          : install dependencies
install:
	OS := $(shell uname)
	ifeq $(OS) Darwin
		sudo gem update —system
		sudo gem install -n /usr/local/bin/ jemoji
		sudo gem install -n /usr/local/bin/ jekyll
		sudo gem install -n /usr/local/bin/ jekyll-feed
		sudo gem install -n /usr/local/bin/ bundler
	else
		gem install jekyll
		gem install jemoji
		gem install jekyll-feed
		gem install bundler
	endif



