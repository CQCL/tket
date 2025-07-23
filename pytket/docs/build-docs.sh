cp -R pytket-docs-theming/_static .
cp pytket-docs-theming/conf.py .

# Get pytket package version to be used in page title
PYTKET_VERSION="$(pip show pytket| grep Version | awk '{print $2}')"

PACKAGE="pytket $PYTKET_VERSION"

# Build the docs setting the html title to show the correct pytket version.
sphinx-build -W -b html . build -D html_title="$PACKAGE API documentation" || exit 1

# Run link checker
sphinx-build -W -b linkcheck . build || exit 1

sphinx-build -W -v -b coverage . build/coverage || exit 1

# Remove copied files. This ensures reusability.
rm -r _static
rm conf.py
