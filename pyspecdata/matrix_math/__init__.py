"""Matrix math utilities.

This module previously imported the :mod:`nnls` implementation eagerly,
which in turn requires a compiled extension (``_nnls``).  The test suite
only exercises the high level :func:`venk_nnls` helper and does not
actually call it.  Importing the compiled extension unconditionally would
therefore raise an :class:`ImportError` during module import on systems
where the extension is not built.  To allow the rest of the package to be
imported without the optional dependency we perform the import lazily
inside :func:`venk_nnls`.
"""

__all__ = ["venk_nnls"]


def venk_nnls(*args, **kwargs):
    """Dispatch to :func:`nnls.venk_nnls` when available.

    The import is deferred so that tests which do not need the compiled
    ``_nnls`` extension can run without it being present.
    """

    from .nnls import venk_nnls as _venk_nnls

    return _venk_nnls(*args, **kwargs)
