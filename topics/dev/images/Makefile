INPUTS := $(wildcard *.plantuml.txt)
OUTPUTS := $(INPUTS:.txt=.svg)

all: plantuml.jar plantuml_options.txt $(OUTPUTS)

$(OUTPUTS): $(INPUTS)
	java -jar plantuml.jar -c plantuml_options.txt -tsvg $(INPUTS)

plantuml.jar:
	wget http://jaist.dl.sourceforge.net/project/plantuml/plantuml.jar
