#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Activate your Python virtual environment (adjust path if necessary)
source ~/base/bin/activate

# Ensure meson, meson-python, and twine are installed
pip install meson meson-python ninja twine

# Create a build directory and configure the project
echo "Configuring the build..."
meson setup builddir

# Compile the project
echo "Building the project..."
meson compile -C builddir

# Install the package locally (optional)
# meson install -C builddir

# Create a wheel for distribution without isolation
echo "Building the wheel..."
python -m build --wheel --no-isolation

# Check the distribution files
echo "Checking the distribution files..."
twine check dist/*

# Upload to PyPI
echo "Uploading to PyPI..."
twine upload dist/*

echo "Done! Your package has been built and uploaded."
