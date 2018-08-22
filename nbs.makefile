CLEAN_NOTEBOOKS := $(wildcard notebooks/*.ipynb.cln)
NOTEBOOKS := $(CLEAN_NOTEBOOKS:.ipynb.cln=.ipynb)

all: $(NOTEBOOKS)

%.ipynb: %.ipynb.cln
	cp $< $@
