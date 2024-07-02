Pytket Docs

# Install python deps
poetry install

# Build html
poetry run make html

# Serve build
npx serve ./build/html
