NOTEBOOKS := $(wildcard notebooks/*.ipynb)
CLEAN_NOTEBOOKS := $(NOTEBOOKS:.ipynb=.ipynb.cln)

all: $(CLEAN_NOTEBOOKS)

%.ipynb.cln: %.ipynb
	notebooks/clear_output_ipynb.py $<
