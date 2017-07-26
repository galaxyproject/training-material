# Settings
JEKYLL=bundle exec jekyll

default: help

serve: ## run a local server
	${JEKYLL} serve
.PHONY: serve

build: ## build files but do not run a server
	${JEKYLL} build
.PHONY: site

check: ## validate HTML
	bundle exec htmlproofer ./_site
.PHONY: check

clean: ## clean up junk files
	@rm -rf _site
	@rm -rf .sass-cache
	@find . -name .DS_Store -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
.PHONY: clean

install: ## install dependencies
	gem install bundler
	bundle install
	bundle update
.PHONY: install

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help