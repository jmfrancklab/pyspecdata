#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define the Python version and corresponding paths
PYTHON_VERSION=cp310-cp310
PYTHON_BIN=/opt/python/$PYTHON_VERSION/bin
MESON_BIN=$PYTHON_BIN/meson

# Run the Docker container and build inside it
docker run --rm -v $(pwd):/io quay.io/pypa/manylinux2014_x86_64 /bin/bash -c "
    # Configure git to recognize the directory as safe
    git config --global safe.directory '*' && \
    # Upgrade pip and install necessary build tools and numpy
    $PYTHON_BIN/python -m pip install --upgrade pip && \
    $PYTHON_BIN/python -m pip install meson meson-python ninja numpy lmfit sympy && \
    # Add the Python bin directory to PATH
    export PATH=\$PATH:$PYTHON_BIN && \
    # Navigate to the mounted directory
    cd /io && \
    # Clean up any previous build configuration
    rm -rf builddir && \
    $MESON_BIN setup builddir && \
    $MESON_BIN compile -C builddir && \
    # Build the wheel and source tar.gz
    $PYTHON_BIN/python -m build --wheel --sdist --no-isolation && \
    # Repair the wheel to ensure manylinux2014 compliance
    auditwheel repair dist/*.whl -w /io/wheelhouse
"

# Copy the source tar.gz to the wheelhouse directory
cp dist/*.tar.gz wheelhouse/

# Check the distribution files
echo "Checking the distribution files..."
twine check wheelhouse/*

# Upload to PyPI
echo "Uploading to PyPI..."
twine upload wheelhouse/*

echo "Done! Your package has been built and uploaded."
