# Settings
JEKYLL=bundle exec jekyll
CHROME=google-chrome-stable
TUTORIALS=$(shell find _site/training-material -name 'tutorial.html' | sed 's/_site\/training-material\///')
SLIDES=$(shell find _site/training-material -name 'slides.html' | sed 's/_site\/training-material\///')
SLIDES+=$(shell find _site/training-material/*/*/slides/* | sed 's/_site\/training-material\///')
SITE_URL=http://localhost:4000/training-material
PDF_DIR=_pdf
REPO=$(shell echo "$${ORIGIN_REPO:-galaxyproject/training-material}")
BRANCH=$(shell echo "$${ORIGIN_BRANCH:-master}")

ifeq ($(shell uname -s),Darwin)
	CHROME=/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome
endif

default: help

serve: ## run a local server
	${JEKYLL} serve -d _site/training-material
.PHONY: serve

detached-serve: clean ## run a local server in detached mode
	${JEKYLL} serve --detach -d _site/training-material
.PHONY: detached-serve

build: clean ## build files but do not run a server
	${JEKYLL} build -d _site/training-material
.PHONY: build

check-html: build ## validate HTML
	bundle exec htmlproofer \
          --assume-extension \
          --http-status-ignore 405,503,999 \
          --url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/" \
          --url-swap "github.com/galaxyproject/training-material/tree/master:github.com/${REPO}/tree/${BRANCH}" \
          --file-ignore "/.*\/files\/.*/","/.*\/node_modules\/.*/" \
          --allow-hash-href \
        ./_site
.PHONY: check-html

check-slides: build  ## check the markdown-formatted links in slides
	find _site -path "**/slides*.html" \
        | xargs -L 1 -I '{}' sh -c \
        "awesome_bot \
           --allow 405 \
           --allow-redirect \
           --white-list localhost,127.0.0.1,fqdn,vimeo.com,drmaa.com \
           --allow-ssl \
           --allow-dupe \
           --skip-save-results \
         -f {}"
.PHONY: check-slides

check-yaml: ## lint yaml files
	find . -path "**/*.yaml" | xargs -L 1 -I '{}' sh -c "yamllint {}"
.PHONY: check-yaml

check: check-yaml check-html check-slides  ## run all checks
.PHONY: check

check-links-gh-pages:  ## validate HTML on gh-pages branch (for daily cron job)
	bundle exec htmlproofer \
          --assume-extension \
          --http-status-ignore 405,503,999 \
          --url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/" \
          --file-ignore "/.*\/files\/.*/" \
          --allow-hash-href \
        .
	find . -path "**/slides*.html" \
        | xargs -L 1 -I '{}' sh -c \
        "awesome_bot \
           --allow 405 \
           --allow-redirect \
           --white-list localhost,127.0.0.1,fqdn,vimeo.com,drmaa.com \
           --allow-ssl \
           --allow-dupe \
           --skip-save-results \
       -f {}"
.PHONY: check-links-gh-pages

clean: ## clean up junk files
	@rm -rf _site
	@rm -rf .sass-cache
	@find . -name .DS_Store -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
.PHONY: clean

install: ## install dependencies
	npm install decktape
	gem install bundler
	bundle install
.PHONY: install

pdf: detached-serve ## generate the PDF of the tutorials and slides
	mkdir -p _pdf
	@for t in $(TUTORIALS); do \
		name="$(PDF_DIR)/$$(echo $$t | tr '/' '-' | sed -e 's/html/pdf/' -e 's/topics-//' -e 's/tutorials-//')"; \
		${CHROME} \
            --headless \
            --disable-gpu \
            --print-to-pdf="$$name" \
            "$(SITE_URL)/$$t?with-answers" \
            2> /dev/null ; \
    done
	@for s in $(SLIDES); do \
		name="$(PDF_DIR)/$$(echo $$s | tr '/' '-' | sed -e 's/html/pdf/' -e 's/topics-//' -e 's/tutorials-//')"; \
		`npm bin`/decktape \
			automatic \
			"$(SITE_URL)/$$s?with-answers" \
			"$$name" \
            2> /dev/null ; \
	done
	pkill -f jekyll
.PHONY: pdf

annotate:
	python bin/add_galaxy_instance_annotations.py
	python bin/add_galaxy_instance_badges.py

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
