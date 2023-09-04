# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	find build/html/ -type f -name "*.html" | xargs sed -i 's/pytket._tket/pytket/g'
	sed -i 's/pytket._tket/pytket/g' build/html/searchindex.js

# for local builds on mac (mac sed is slightly different)
html-mac:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	find build/html/ -type f -name "*.html" | xargs sed -e 's/pytket._tket/pytket/g' -i ""
	sed -e 's/pytket._tket/pytket/g' -i "" build/html/searchindex.js
