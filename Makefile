# Settings
JEKYLL=jekyll
PORT?=4000
HOST?=0.0.0.0
FLAGS?=""
ENV?="development"
CHROME=google-chrome-stable
PDF_HOST?=127.0.0.1
SITE_URL=http://${PDF_HOST}:${PORT}/training-material
PDF_DIR=_pdf
REPO=$(shell echo "$${ORIGIN_REPO:-galaxyproject/training-material}")
BRANCH=$(shell echo "$${ORIGIN_BRANCH:-master}")
MINICONDA_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
SHELL=bash
RUBY_VERSION=2.4.4
CONDA_ENV=galaxy_training_material

ifeq ($(shell uname -s),Darwin)
	CHROME=/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome
	MINICONDA_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
endif

CONDA=$(shell which conda)
ifeq ($(CONDA),)
	CONDA=${HOME}/miniconda3/bin/conda
endif

default: help

install-conda: ## install Miniconda
	curl -L $(MINICONDA_URL) -o miniconda.sh
	bash miniconda.sh -b
.PHONY: install-conda

create-env: ## create conda environment
	if ${CONDA} env list | grep '^${CONDA_ENV}'; then \
	    ${CONDA} env update -f environment.yml; \
	else \
	    ${CONDA} env create -f environment.yml; \
	fi
.PHONY: create-env

ACTIVATE_ENV = source $(shell dirname $(dir $(CONDA)))/bin/activate $(CONDA_ENV)

install: clean create-env ## install dependencies
	$(ACTIVATE_ENV) && \
		gem update --system && \
		gem install addressable:'2.5.2' jekyll jekyll-feed jekyll-scholar jekyll-redirect-from jekyll-last-modified-at csl-styles awesome_bot html-proofer pkg-config kwalify
.PHONY: install

serve: ## run a local server (You can specify PORT=, HOST=, and FLAGS= to set the port, host or to pass additional flags)
	@echo "Tip: Want faster builds? Use 'serve-quick' in place of 'serve'."
	@echo "Tip: to serve in incremental mode (faster rebuilds), use the command: make serve FLAGS=--incremental" && echo "" && \
	$(ACTIVATE_ENV) && \
		mv Gemfile Gemfile.backup || true && \
		mv Gemfile.lock Gemfile.lock.backup || true && \
		${JEKYLL} serve --strict_front_matter -d _site/training-material -P ${PORT} -H ${HOST} ${FLAGS}
.PHONY: serve

serve-quick: ## run a local server (faster, some plugins disabled for speed)
	@echo "This will build the website with citations and other content disabled, and incremental on by default. To run the full preview (slower), use make serve" && echo "" && \
	$(ACTIVATE_ENV) && \
		mv Gemfile Gemfile.backup || true && \
		mv Gemfile.lock Gemfile.lock.backup || true && \
		${JEKYLL} serve --strict_front_matter -d _site/training-material --incremental --config _config.yml,_config-dev.yml -P ${PORT} -H ${HOST} ${FLAGS}
.PHONY: serve-quick

build: clean ## build files but do not run a server (You can specify FLAGS= to pass additional flags to Jekyll)
	$(ACTIVATE_ENV) && \
		mv Gemfile Gemfile.backup || true && \
		mv Gemfile.lock Gemfile.lock.backup || true && \
		JEKYLL_ENV=${ENV} ${JEKYLL} build --strict_front_matter -d _site/training-material ${FLAGS}
.PHONY: build

check-frontmatter: ## Validate the frontmatter
	$(ACTIVATE_ENV) && \
		find topics/ -name tutorial.md -or -name slides.html -or -name metadata.yaml | \
	    xargs -n1 ruby bin/validate-frontmatter.rb
.PHONY: check-frontmatter

check-html: build ## validate HTML
	$(ACTIVATE_ENV) && \
	  	htmlproofer \
	      	--assume-extension \
	      	--http-status-ignore 405,503,999 \
	      	--url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/" \
	      	--url-swap "github.com/galaxyproject/training-material/tree/master:github.com/${REPO}/tree/${BRANCH}" \
	      	--file-ignore "/.*\/files\/.*/","/.*\/node_modules\/.*/" \
	      	--allow-hash-href \
	      	./_site
.PHONY: check-html

check-workflows: ## validate Workflows
	$(ACTIVATE_ENV) && \
		bash bin/validate-json.sh && \
		bash bin/validate-workflow-tags.sh
.PHONY: check-workflows

check-references: build ## validate no missing references
	$(ACTIVATE_ENV) && \
		bash bin/validate-references.sh
.PHONY: check-references

check-html-internal: build ## validate HTML (internal links only)
	$(ACTIVATE_ENV) && \
		htmlproofer \
	      	--assume-extension \
	      	--http-status-ignore 405,503,999 \
	      	--url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/","/.*slides.html#/" \
	      	--url-swap "github.com/galaxyproject/training-material/tree/master:github.com/${REPO}/tree/${BRANCH}" \
	      	--file-ignore "/.*\/files\/.*/","/.*\/node_modules\/.*/" \
	      	--disable-external \
	      	--allow-hash-href \
	      	./_site
.PHONY: check-html-internal

check-slides: build  ## check the markdown-formatted links in slides
	$(ACTIVATE_ENV) && \
		find _site -path "**/slides*.html" \
	      	| xargs -L 1 -I '{}' sh -c "awesome_bot \
				--allow 405 \
				--allow-redirect \
				--white-list localhost,127.0.0.1,fqdn,vimeo.com,drmaa.com \
				--allow-ssl \
				--allow-dupe \
				--skip-save-results \
				-f {}"
.PHONY: check-slides

check-yaml: ## lint yaml files
	$(ACTIVATE_ENV) && \
		find . -name "*.yaml" | grep -v .github | xargs -L 1 -I '{}' sh -c "yamllint {}" \
		find topics -name '*.yml' | xargs -L 1 -I '{}' sh -c "yamllint {}" \
		ruby bin/check-contributors.rb
.PHONY: check-yaml

check-snippets: ## lint snippets
	./bin/check-for-trailing-newline
.PHONY: check-snippets

check-framework:
	$(ACTIVATE_ENV) && \
		ruby _plugins/jekyll-notranslate.rb
.PHONY: check-framework

check-broken-boxes: build ## List tutorials containing broken boxes
	./bin/check-broken-boxes
.PHONY: check-broken-boxes

check: check-yaml check-frontmatter check-html-internal check-html check-broken-boxes check-slides check-workflows check-references check-snippets ## run all checks
.PHONY: check

lint: ## run all linting checks
	$(MAKE) check-yaml
	$(MAKE) check-frontmatter
	$(MAKE) check-workflows
	$(MAKE) check-references
	$(MAKE) check-snippets
.PHONY: lint

check-links-gh-pages:  ## validate HTML on gh-pages branch (for daily cron job)
	$(ACTIVATE_ENV) && \
	  	htmlproofer \
			--assume-extension \
			--http-status-ignore 405,503,999 \
			--url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/" \
			--file-ignore "/.*\/files\/.*/" \
			--allow-hash-href \
			. && \
		find . -path "**/slides*.html" \
			| xargs -L 1 -I '{}' sh -c "awesome_bot \
				--allow 405 \
				--allow-redirect \
				--white-list localhost,127.0.0.1,fqdn,vimeo.com,drmaa.com \
				--allow-ssl \
				--allow-dupe \
				--skip-save-results \
				-f {}"
.PHONY: check-links-gh-pages

TUTORIAL_PDFS=$(shell find _site/training-material -name 'tutorial.html' | sed 's/html$$/pdf/g')
SLIDE_PDFS=$(shell find _site/training-material -name 'slides.html' | sed 's/html$$/pdf/g')
SLIDE_PDFS+=$(shell find _site/training-material/*/*/slides/* | sed 's/html$$/pdf/g')

pdf: $(SLIDE_PDFS) $(TUTORIAL_PDFS) ## generate the PDF of the tutorials and slides
.PHONY: pdf

_site/%/tutorial.pdf: _site/%/tutorial.html
	if ! grep 'http-equiv="refresh"' $< --quiet; then \
		$(ACTIVATE_ENV) && \
		sed "s|/training-material/|$(shell pwd)/_site/training-material/|g" $< | \
		sed "s|<head>|<head><base href=\"file://$(shell pwd)/$(<:_site/training/material%=%)\">|" | \
		wkhtmltopdf \
		    --enable-javascript --javascript-delay 1000 \
			- $@; \
	fi

_site/%.pdf: _site/%.html
	if ! grep 'http-equiv="refresh"' $< --quiet; then \
		$(ACTIVATE_ENV) && \
		sed "s|/training-material/|$(shell pwd)/_site/training-material/|g" $< | \
		sed "s|<head>|<head><base href=\"file://$(shell pwd)/$(<:_site/training/material%=%)\">|" | \
		wkhtmltopdf \
		    --enable-javascript --javascript-delay 3000 --page-width 700px --page-height 530px -B 5px -L 5px -R 5px -T 5px \
			--user-style-sheet bin/slides-fix.css \
			- $@; \
	fi

AWS_UPLOAD?=""
VIDEOS := $(shell find topics -name 'slides.html' | xargs ./bin/filter-has-videos)
video: $(VIDEOS:topics/%.html=_site/training-material/topics/%.mp4) ## Build videos where possible

_site/training-material/%/slides.mp4: _site/training-material/%/slides.pdf %/slides.html
	./bin/ari.sh $^ $@ $(AWS_UPLOAD)

annotate: ## annotate the tutorials with usable Galaxy instances and generate badges
	${ACTIVATE_ENV} && \
	bash bin/workflow_to_tool_yaml.sh && \
	python bin/add_galaxy_instance_annotations.py && \
	python bin/add_galaxy_instance_badges.py
.PHONY: annotate

clean: ## clean up junk files
	@rm -rf _site
	@rm -rf .sass-cache
	@rm -rf .bundle
	@rm -rf vendor
	@rm -rf node_modules
	@rm -rf .jekyll-cache
	@rm -rf .jekyll-metadata
	@find . -name .DS_Store -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
.PHONY: clean

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
