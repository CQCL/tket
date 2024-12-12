cp -R pytket-docs-theming/_static .
cp -R pytket-docs-theming/quantinuum-sphinx .
cp pytket-docs-theming/conf.py .

sphinx-build -b html . build -D html_title="pytket"
find build/ -type f -name "*.html" | xargs sed -e 's/pytket._tket/pytket/g' -i ""
sed -i '' 's/pytket._tket/pytket/g' build/searchindex.js

# Remove copied files. This ensures reusability.
rm -r _static 
rm -r quantinuum-sphinx
rm conf.py