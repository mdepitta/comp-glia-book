# -*- Mode: Makefile; indent-tabs-mode: t; tab-width: 8 -*-

.PHONY: all

# Code examples/snippets
code_files = example_1_COBA.py example_2_gchi_astrocyte.py example_3_io_synapse.py \
             example_4_synrel.py example_4_rsmean.py example_5_astro_ring.py \
             example_6_COBA_with_astro.py
base_code_files = $(foreach file,$(code_files),code/$(file))
clean_code_files = $(foreach file,$(code_files),code/clean_examples/$(file))

# Example 1 (COBA)
example_1_fig_files = text/figures/results/example_1_COBA_Figure_1.eps \
                      text/figures/results/example_1_COBA_Figure_2.eps
example_2_fig_files = text/figures/results/example_2_gchi_astrocyte_Figure.eps
example_3_fig_files = text/figures/results/example_3_io_synapse_Figure.eps
example_4a_fig_files = text/figures/results/example_4_synrel_Figure.eps
example_4b_fig_files = text/figures/results/example_4_rsmean_Figure.eps
example_5_fig_files = text/figures/results/example_5_astro_ring_Figure.eps
example_6_fig_files = text/figures/results/example_6_COBA_with_astro_Figure.eps

example_figures = $(example_1_fig_files) $(example_2_fig_files) $(example_3_fig_files) \
                  $(example_4a_fig_files) $(example_4b_fig_files) $(example_5_fig_files) \
                  $(example_6_fig_files)

all : text/chapter18_brian.pdf text/parameters.pdf

# Main document
text/chapter18_brian.pdf : text/chapter18_brian.tex text/references.bib $(base_code_files) $(clean_code_files) $(example_figures)
	cd text; pdflatex chapter18_brian.tex
	cd text; bibtex chapter18_brian
	cd text; pdflatex chapter18_brian.tex && pdflatex chapter18_brian.tex

# Parameter overview PDF
text/parameters.pdf : text/parameters.tex
	cd text; pdflatex parameters.tex

text/parameters.tex : $(base_code_files) code/extract_parameters.py
	cd code; python extract_parameters.py

# Code examples/snippets
$(clean_code_files) : $(base_code_files)
	cd code; python extract_code.py

# Figures (i.e., results of example runs)
$(example_1_fig_files) : code/example_1_COBA.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_1_COBA.py

$(example_2_fig_files) : code/example_2_gchi_astrocyte.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_2_gchi_astrocyte.py

$(example_3_fig_files) : code/example_3_io_synapse.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_3_io_synapse.py

$(example_4a_fig_files) : code/example_4_synrel.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_4_synrel.py

$(example_4b_fig_files) : code/example_4_rsmean.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_4_rsmean.py

$(example_5_fig_files) : code/example_5_astro_ring.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_5_astro_ring.py

$(example_6_fig_files) : code/example_6_COBA_with_astro.py code/figures.mplstyle code/plot_utils.py
	cd code; python example_6_COBA_with_astro.py

