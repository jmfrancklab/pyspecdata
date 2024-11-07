#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define the Python version and corresponding paths
PYTHON_VERSION=cp311-cp311
PYTHON_BIN=/opt/python/$PYTHON_VERSION/bin
MESON_BIN=$PYTHON_BIN/meson

# Run the Docker container and build inside it
docker run --rm -v $(pwd):/io quay.io/pypa/manylinux2014_x86_64 /bin/bash -c "
    # Upgrade pip and install necessary build tools
    $PYTHON_BIN/python -m pip install --upgrade pip && \
    $PYTHON_BIN/python -m pip install meson meson-python ninja && \
    # Add the Python bin directory to PATH
    export PATH=\$PATH:$PYTHON_BIN && \
    # Navigate to the mounted directory
    cd /io && \
    # Clean up any previous build configuration
    rm -rf builddir && \
    $MESON_BIN setup builddir && \
    $MESON_BIN compile -C builddir && \
    # Build the wheel
    $PYTHON_BIN/python -m build --wheel --no-isolation && \
    # Repair the wheel to ensure manylinux2014 compliance
    auditwheel repair dist/*.whl -w /io/wheelhouse
"

# Check the distribution files
echo "Checking the distribution files..."
twine check wheelhouse/*

# Upload to PyPI
echo "Uploading to PyPI..."
twine upload wheelhouse/*

echo "Done! Your package has been built and uploaded."
