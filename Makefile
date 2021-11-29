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
BRANCH=$(shell echo "$${ORIGIN_BRANCH:-main}")
MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
SHELL=bash
RUBY_VERSION=2.4.4
CONDA_ENV=galaxy_training_material

ifeq ($(shell uname -s),Darwin)
	CHROME=/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome
	MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
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
		gem update --no-document --system && \
		ICONV_LIBS="-L${CONDA_PREFIX}/lib/ -liconv" gem install --no-document addressable:'2.5.2' jekyll jekyll-feed jekyll-redirect-from jekyll-last-modified-at csl-styles awesome_bot html-proofer pkg-config kwalify bibtex-ruby citeproc-ruby
.PHONY: install

bundle-install: clean  ## install gems if Ruby is already present (e.g. on gitpod.io)
	bundle install
.PHONE: bundle-install

serve: api/swagger.json ## run a local server (You can specify PORT=, HOST=, and FLAGS= to set the port, host or to pass additional flags)
	@echo "Tip: Want faster builds? Use 'serve-quick' in place of 'serve'."
	@echo "Tip: to serve in incremental mode (faster rebuilds), use the command: make serve FLAGS=--incremental" && echo "" && \
	$(ACTIVATE_ENV) && \
		mv Gemfile Gemfile.backup || true && \
		mv Gemfile.lock Gemfile.lock.backup || true && \
		${JEKYLL} serve --strict_front_matter -d _site/training-material -P ${PORT} -H ${HOST} ${FLAGS}
.PHONY: serve

serve-quick: api/swagger.json ## run a local server (faster, some plugins disabled for speed)
	@echo "This will build the website with citations and other content disabled, and incremental on by default. To run the full preview (slower), use make serve" && echo "" && \
	$(ACTIVATE_ENV) && \
		mv Gemfile Gemfile.backup || true && \
		mv Gemfile.lock Gemfile.lock.backup || true && \
		${JEKYLL} serve --strict_front_matter -d _site/training-material --incremental --config _config.yml,_config-dev.yml -P ${PORT} -H ${HOST} ${FLAGS}
.PHONY: serve-quick

serve-gitpod: bundle-install api/swagger.json  ## run a server on a gitpod.io environment
	bundle exec jekyll serve --config _config.yml --incremental
.PHONY: serve-gitpod

build-gitpod: bundle-install api/swagger.json  ## run a build on a gitpod.io environment
	bundle exec jekyll build --config _config.yml
.PHONY: build-gitpod

build: clean api/swagger.json ## build files but do not run a server (You can specify FLAGS= to pass additional flags to Jekyll)
	$(ACTIVATE_ENV) && \
		mv Gemfile Gemfile.backup || true && \
		mv Gemfile.lock Gemfile.lock.backup || true && \
		JEKYLL_ENV=${ENV} ${JEKYLL} build --strict_front_matter -d _site/training-material ${FLAGS}
.PHONY: build

check-frontmatter: ## Validate the frontmatter
	$(ACTIVATE_ENV) && \
		bundle exec ruby bin/validate-frontmatter.rb
.PHONY: check-frontmatter

check-contributors: ## Validate the contributors.yaml file
	$(ACTIVATE_ENV) && \
		bundle exec ruby bin/validate-contributors.rb
.PHONY: check-contributors

_check-html: # Internal
	$(ACTIVATE_ENV) && \
	  	htmlproofer \
	      	--assume-extension \
	      	--http-status-ignore 405,503,999 \
	      	--url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/" \
	      	--url-swap "github.com/galaxyproject/training-material/tree/main:github.com/${REPO}/tree/${BRANCH}" \
	      	--file-ignore "/.*\/files\/.*/","/.*\/node_modules\/.*/","/\/tutorials\/.*\/docker\//" \
	      	--allow-hash-href \
	      	./_site
.PHONY: _check-html

check-html: build ## validate HTML
	$(MAKE) _check-html
.PHONY: check-html

check-workflows: ## validate Workflows
	find topics -name '*.ga' | grep /workflows/ | xargs -P8 -n1 bash bin/validate-workflow.sh
.PHONY: check-workflows

check-references: build ## validate no missing references
	bash bin/validate-references.sh
.PHONY: check-references

_check-html-internal: # Internal
	$(ACTIVATE_ENV) && \
		htmlproofer \
	      	--assume-extension \
	      	--http-status-ignore 405,503,999 \
	      	--url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/","/.*slides.html#/","/#embedded_jbrowse/","/.*videos.*.mp4.png/" \
	      	--url-swap "github.com/galaxyproject/training-material/tree/main:github.com/${REPO}/tree/${BRANCH}" \
	      	--file-ignore "/.*\/files\/.*/","/.*\/node_modules\/.*/","/\/tutorials\/.*\/docker\//" \
	      	--disable-external \
	      	--allow-hash-href \
	      	./_site
.PHONY: _check-html-internal

check-html-internal: build ## validate HTML (internal links only)
	$(MAKE) _check-html-internal
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
	find . -name '*.yaml' | grep -v .github | xargs -L 1 -I '{}' sh -c "yamllint -c .yamllint {}"
.PHONY: check-yaml

check-diffs: ## lint diffs in tutorials
	find ./topics/admin/ -name '*.md' -type f -print0 | xargs -0 python bin/lint-diffs.py
.PHONY: check-diffs

check-tool-links: ## lint tool links
	@bash ./bin/check-broken-tool-links.sh
.PHONY: check-tool-links

check-framework:
	$(ACTIVATE_ENV) && \
		ruby _plugins/jekyll-notranslate.rb
.PHONY: check-framework

check-broken-boxes: build ## List tutorials containing broken boxes
	./bin/check-broken-boxes
.PHONY: check-broken-boxes

check: check-html-internal check-html check-broken-boxes check-slides ## run checks which require compiled HTML
.PHONY: check

lint: check-frontmatter check-workflows check-tool-links ## run linting checks which do not require a built site
.PHONY: lint

check-links-gh-pages:  ## validate HTML on gh-pages branch (for daily cron job)
	$(ACTIVATE_ENV) && \
	  	htmlproofer \
			--assume-extension \
			--http-status-ignore 405,503,999 \
			--url-ignore "/.*localhost.*/","/.*vimeo\.com.*/","/.*gitter\.im.*/","/.*drmaa\.org.*/" \
			--file-ignore "/.*\/files\/.*/","/\/tutorials\/.*\/docker\//" \
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


_site/%/introduction.pdf: _site/%/introduction.html
	$(ACTIVATE_ENV) && \
	$(shell npm bin)/http-server _site -p 9876 & \
	$(shell npm bin)/decktape automatic -s 1920x1080 http://localhost:9876/$(<:_site/%=%) $@; \

_site/%/slides.pdf: _site/%/slides.html
	$(ACTIVATE_ENV) && \
	$(shell npm bin)/http-server _site -p 9876 & \
	$(shell npm bin)/decktape automatic -s 1920x1080 http://localhost:9876/$(<:_site/%=%) $@; \

_site/%/slides_ES.pdf: _site/%/slides_ES.html
	$(ACTIVATE_ENV) && \
	$(shell npm bin)/http-server _site -p 9876 & \
	$(shell npm bin)/decktape automatic -s 1920x1080 http://localhost:9876/$(<:_site/%=%) $@; \

_site/%/slides_CAT_ES.pdf: _site/%/slides_CAT_ES.html
	$(ACTIVATE_ENV) && \
	$(shell npm bin)/http-server _site -p 9876 & \
	$(shell npm bin)/decktape automatic -s 1920x1080 http://localhost:9876/$(<:_site/%=%) $@; \

video: ## Build all videos
	bash bin/ari-make.sh

annotate: ## annotate the tutorials with usable Galaxy instances and generate badges
	${ACTIVATE_ENV} && \
	bash bin/workflow_to_tool_yaml.sh && \
	python bin/add_galaxy_instance_annotations.py && \
	python bin/add_galaxy_instance_badges.py
.PHONY: annotate

rebuild-search-index: ## Rebuild search index
	node bin/lunr-index.js > search.json

api/swagger.json: metadata/swagger.yaml
	cat metadata/swagger.yaml | python bin/yaml2json.py > api/swagger.json

clean: ## clean up junk files
	@rm -rf _site
	@rm -rf .sass-cache
	@rm -rf .bundle
	@rm -rf vendor
	@rm -rf node_modules
	@rm -rf .jekyll-metadata
	@find . -name .DS_Store -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
.PHONY: clean

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
