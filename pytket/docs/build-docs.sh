cp -R pytket-docs-theming/_static .
cp -R pytket-docs-theming/quantinuum-sphinx .
cp pytket-docs-theming/conf.py .

# Get pytket package version to be used in page title
PYTKET_VERSION="$(pip show pytket| grep Version | awk '{print $2}')"

PACKAGE="pytket $PYTKET_VERSION"

# Build the docs setting the html title to show the correct pytket version.
sphinx-build -b html . build -D html_title="$PACKAGE API documentation"

sphinx-build -b coverage . build/coverage

# Replace unnecessary _tket for all classes and functions in the built html
# Apple MACOSX and Linux have differing sed replace syntax
if [[ "$OSTYPE" == "darwin"* ]]; then
    find build/ -type f -name "*.html" | xargs sed -e 's/pytket._tket/pytket/g' -i ""
    sed -i '' 's/pytket._tket/pytket/g' build/searchindex.js
else
    find build/ -type f -name "*.html" | xargs sed -i 's/pytket._tket/pytket/g'
    sed -i 's/pytket._tket/pytket/g' build/searchindex.js
fi

# Remove copied files. This ensures reusability.
rm -r _static 
rm -r quantinuum-sphinx
rm conf.py