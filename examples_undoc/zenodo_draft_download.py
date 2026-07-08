"""Download a file from an unpublished Zenodo draft deposition
---------------------------------------------------------------

Published Zenodo record downloads do not require authentication.  Draft
deposition downloads select a named token from ``[zenodo_tokens]`` in
``~/.pyspecdata``::

    [zenodo]
    default_token = lab

    [zenodo_tokens]
    lab = /path/to/zenodo.token

Each token file should contain a Zenodo personal access token with appropriate
deposit permissions.  Keep token files local and never commit them to git.
Supplying ``token_name`` tells :func:`zenodo_download` to use the authenticated
draft-deposition API.  See https://developers.zenodo.org/#rest-api for API and
authentication details.
"""

from pyspecdata import zenodo_download


path = zenodo_download(
    "123456",
    "my_file.h5",
    exp_type="ODNP",
    token_name="lab",
)
print(f"Downloaded to {path}")
