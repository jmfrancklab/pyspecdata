"""Download a file from Zenodo using :func:`search_filename`
-----------------------------------------------------------

This example shows how to retrieve a file from Zenodo when it is not already
present locally.  The download occurs transparently via the ``zenodo``
keyword of :func:`pyspecdata.search_filename`.
"""

from pyspecdata import search_filename

path = search_filename(
    "testdata.h5",
    exp_type="example_data",
    unique=True,
    zenodo="https://zenodo.org/record/123456/files/testdata.h5?download=1",
)

print(f"Downloaded to {path}")
