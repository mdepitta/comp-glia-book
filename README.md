# Overview

The repo contains code used for simulations and figures in the book "Computational Glioscience," De Pitta' M. and Berry, H. eds., Springer (expected 2018). The repo contains different folders accoridng to the chapters in the book. Content of each folder is detailed as following.

#--------------------------------------------------------
# Ch19.Stimberg
#--------------------------------------------------------

The code directory contains one file for each example that can be run independently and by default creates and saves images to the text/figures/results directory. It also contains a script named "extract_code.py" that takes the code files in the code directory (which are the only ones that should ever be edited) and creates "clean examples" from them (that we'll can later include in Brian or put on a website or something), as well as the code snippets that will be included in the text. See below for notes on how this script works.

The LaTeX file needs that extract_code.py has been run and that all the figures have been generated. The easiest way to make sure that everything is up-to-date is to run "make" in the main directory, it will re-run everything that is missing, build the PDF, call bibtex, etc.

The code contained in this folder builds the whole chapter from the mere TEX file, running all simulations and producing clean Brian2 examples. Once downloaded, on the command line in the local repo folder, type 'cd Ch19.Stimberg | ./make'.

## Code annotations
In the code, there are a few strange looking comments. These are used so that we can use them as a source for the code snippets included in the paper without having to worry that the two codes get out-of-sync, i.e. that the code used to generate the figures does not correspond to the code we are presenting in the paper. There are a few comments that are interpreted specially:

* DELETE: This will completely delete this line when creating the "clean examples". This should be for stuff that is paper-specific, e.g. where to save the figure files too and maybe some very specific customisations of the plot
* INCLUDE BEGIN and INCLUDE END: This will start/end a block that will be used as a code snippet in the paper. Blocks will be numbered starting at 1, so the first such block of `example_1_COBA.py` can be included with `\lstinputlisting{code/example_1_COBA_1.py}` in the LaTeX document
* ELLIPSIS BEGIN and ELLIPSIS END: This will replace the respective block (that has to be within an include block) by "# ..." in the generated code snippet. This is useful when we show some code (e.g. astrocyte equations) that has been shown before but has some new additions.
