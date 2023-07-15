#!/usr/bin/bash

# build and open docs locally
# (only works on windows -- on Mac, I think you say "open" instead of "start")

# you need these packages
#
# conda -f conda-forge make
# pip install sphinx-gallery

if grep -q '\bquit\b' -r ../examples/
then
    echo "example files may not use the quit statement!!!!"
    echo ""
    echo "the following files use quit:"
    grep -e '\bquit\b' -rl ../examples/
else
    make html
    start _build/html/auto_examples/index.html
fi
