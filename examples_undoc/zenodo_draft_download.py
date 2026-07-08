"""Download a file from an unpublished Zenodo draft deposition.

Published Zenodo record downloads do not require authentication.  Draft
deposition downloads reuse the existing ``[zenodo] token_file`` setting.  The
token file should contain a Zenodo personal access token with appropriate
deposit permissions.  Keep that file local and never commit it to git.
"""

from pyspecdata import zenodo_download


path = zenodo_download(
    "123456",
    "my_file.h5",
    exp_type="ODNP",
    draft=True,
)
print(f"Downloaded to {path}")
