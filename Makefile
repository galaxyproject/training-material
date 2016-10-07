DECKTAPE_DIR := decktape-1.0.0/
OS := $(shell uname)
SLIDE_PATTERN := slides/

pdf:
    ifeq ($(URL), )
		@echo "Please specify a file. 'make pdf FILE=...'"
    else
		$(DECKTAPE_DIR)phantomjs $(DECKTAPE_DIR)decktape.js $(FILE) $(basename $(FILE)).pdf
    endif

pdf-recursive:
    ifeq ($(DIR), )
		@echo "Please specify a directory. 'make pdf-recursice DIR=...'"
    else
		$(eval FILES := $(shell find $(DIR) -name *.html -print | grep -F '$(SLIDE_PATTERN)'))
		$(foreach var, $(FILES), $(shell $(DECKTAPE_DIR)phantomjs \
			$(DECKTAPE_DIR)decktape.js $(var) $(basename $(var)).pdf >&2))
    endif

install-decktape:
	curl -L https://github.com/astefanutti/decktape/archive/v1.0.0.tar.gz | tar -xz --exclude phantomjs
    ifeq ($(OS),Darwin)
		curl -L https://github.com/astefanutti/decktape/releases/download/v1.0.0/phantomjs-osx-cocoa-x86-64 \
		-o $(DECKTAPE_DIR)phantomjs
    else
		curl -L https://github.com/astefanutti/decktape/releases/download/v1.0.0/phantomjs-linux-x86-64 \
		-o $(DECKTAPE_DIR)phantomjs
    endif
	chmod +x $(DECKTAPE_DIR)phantomjs
