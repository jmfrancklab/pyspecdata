import logging
import numpy as np
import time
from numpy import newaxis
from ..general_functions import strm

logger = logging.getLogger("pyspecdata.matrix_math")


def along(self, dimname, rename_redundant=None):
    r"""Specifies the dimension for the next matrix
    multiplication (represents the rows/columns).

    Parameters
    ==========
    dimname: str
        The next time matrix multiplication is called,
        'dimname' will be summed over.
        That is, `dimname` will become the columns position if this
        is the first matrix.

        If `along` is not called for the second matrix,
        `dimname` will also take the position of rows for that
        matrix!
    rename_redundant: tuple of str or (default) None
        If you are multiplying two different matrices,
        then it is only sensible that before the multiplication,
        you should identify the dimension representing the row
        space of the right matrix and the column space of the left
        matrix with different names.

        **However** sometimes
        (*e.g.* constructing projection matrices)
        you may want to start with two matrices where both the
        row space of the right matrix and the column space of the
        left have the same name.
        If so, you will want to rename the column space of the
        resulting matrix -- then you pass
        ``rename_redundant=('orig name','new name')``
    """
    assert dimname in self.dimlabels
    self._matmul_along = dimname
    if rename_redundant is None:
        if hasattr(self, "_rename_redundant"):
            del self._rename_redundant
    else:
        self._rename_redundant = rename_redundant
        assert self._rename_redundant[0] in self.dimlabels
    return self
def dot(self, arg):
    """This will perform a dot product or a matrix multiplication.
    If one dimension in ``arg`` matches that in ``self``,
    it will dot along that dimension
    (take a matrix multiplication where that
    dimension represents the columns of
    ``self`` and the rows of ``arg``)

    Note that if you have your dimensions
    named "rows" and "columns", this will be
    very confusing, but if you have your
    dimensions named in terms of the vector
    basis they are defined/live in, this
    makes sense.

    If there are zero or no matching
    dimensions, then use
    :func:`~pyspecdata.nddata.along` to
    specify the dimensions for matrix
    multiplication / dot product.

    .. literalinclude:: ../examples/matrix_mult.py
    """
    time_dotstart = time.time()
    if hasattr(self, "_matmul_along"):
        dot_dim_self = self._matmul_along
        if hasattr(arg, "_matmul_along"):
            dot_dim_arg = arg._matmul_along
        else:
            dot_dim_arg = dot_dim_self
    elif hasattr(arg, "_matmul_along"):
        dot_dim_arg = arg._matmul_along
        dot_dim_self = dot_dim_arg
    else:
        matching_dims = set(self.dimlabels) & set(arg.dimlabels)
        if len(matching_dims) == 1:
            dot_dim_self = list(matching_dims)[0]
            dot_dim_arg = dot_dim_self
        else:
            raise ValueError(
                "I can't determine the dimension for matrix"
                " multiplication!  Either have only one matching dimension,"
                " or use .along("
            )
    logger.debug(
        strm(
            "initial shapes are self:",
            self.shape,
            "arg:",
            arg.shape,
            "dotdims (respectively)",
            dot_dim_self,
            dot_dim_arg,
        )
    )
    # {{{ unset the "along" setting
    if hasattr(self, "_matmul_along"):
        del self._matmul_along
    if hasattr(arg, "_matmul_along"):
        del arg._matmul_along
    # }}}
    if hasattr(self, "_rename_redundant"):
        rename_redundant = self._rename_redundant
        del self._rename_redundant
    else:
        rename_redundant = ["XXXRENAMED", "XXXRENAMED"]
    uninvolved_dims = [
        [
            rename_redundant[1] if j == rename_redundant[0] else j
            for j in self.dimlabels
            if j != dot_dim_self
        ]
    ]
    uninvolved_dims += [[j for j in arg.dimlabels if j != dot_dim_arg]]
    mult_dims = [[dot_dim_self], [dot_dim_arg]]
    # {{{ dimensions involved in multiplication must not be
    # repeated
    for k, l in [(1, 0), (0, 1)]:
        for j in range(len(uninvolved_dims[k]) - 1, -1, -1):
            if uninvolved_dims[k][j] not in uninvolved_dims[l]:
                mult_dims[k].insert(k, uninvolved_dims[k].pop(j))
                break
    orig_mult_dims = mult_dims[0] + mult_dims[1]
    logger.debug(strm("mult_dims", mult_dims))
    for j in range(2):
        if len(mult_dims[j]) < 2:
            mult_dims[j].insert(j, "XXX_ADDED_XXX")
    logger.debug(strm("mult_dims", mult_dims))
    # }}}
    # {{{ construct the dictionary right away
    axis_coords_dict = self.mkd(self.axis_coords)
    axis_units_dict = self.mkd(self.axis_coords_units)
    axis_coords_error_dict = self.mkd(self.axis_coords_error)
    arg_axis_coords_dict = arg.mkd(arg.axis_coords)
    arg_axis_units_dict = arg.mkd(arg.axis_coords_units)
    arg_axis_coords_error_dict = arg.mkd(arg.axis_coords_error)
    shared_info = set(self.dimlabels) & set(arg.dimlabels)
    for j in shared_info:
        if axis_coords_dict[j] is None:
            assert arg_axis_coords_dict[j] is None
        else:
            assert all(axis_coords_dict[j] == arg_axis_coords_dict[j])
        assert axis_units_dict[j] == arg_axis_units_dict[j]
        if axis_coords_error_dict[j] is None:
            assert arg_axis_coords_error_dict[j] is None
        else:
            assert all(axis_coords_error_dict[j] == arg_axis_coords_error_dict[j])
    info_needed_from_arg = (
        (
            set(uninvolved_dims[0])
            | set(uninvolved_dims[1])
            | set(orig_mult_dims)
        )
        - set(self.dimlabels)
        - set([rename_redundant[1]])
    )
    logger.debug(
        strm("info needed from arg for dimensions", info_needed_from_arg)
    )
    for j in info_needed_from_arg:
        axis_coords_dict[j] = arg_axis_coords_dict[j]
        axis_units_dict[j] = arg_axis_units_dict[j]
        axis_coords_error_dict[j] = arg_axis_coords_error_dict[j]
    # }}}
    if (self.get_error() is not None) or (arg.get_error() is not None):
        raise ValueError(
            "we plan to include error propagation here, but not yet provided"
        )
    # we need to get into the shape:
    # uninvolved_dims+mult_dims
    output_shape = [mult_dims[0][0], mult_dims[1][1]]
    # {{{ pull alternately from the end of arg and self, matching
    # where possible
    self_temp_list = list(uninvolved_dims[0])
    arg_temp_list = list(uninvolved_dims[1])
    while len(self_temp_list) > 0 or len(arg_temp_list) > 0:
        if len(arg_temp_list) > 0:
            thisdim = arg_temp_list.pop(-1)
            if thisdim in self_temp_list:
                self_temp_list.remove(thisdim)
            output_shape = [thisdim] + output_shape
        if len(self_temp_list) > 0:
            thisdim = self_temp_list.pop(-1)
            if thisdim in arg_temp_list:
                arg_temp_list.remove(thisdim)
            output_shape = [thisdim] + output_shape
    # }}}
    logger.debug(strm("output_shape", output_shape))
    self_new_order = tuple(
        self.axn(rename_redundant[0] if j == rename_redundant[1] else j)
        for j in output_shape[:-1] + [dot_dim_self]
        if j in self.dimlabels + [rename_redundant[1]]
    )
    logger.debug(
        strm(
            "self_new_order",
            self_new_order,
            "self.data.shape",
            self.data.shape,
        )
    )
    self_formult_data = self.data.transpose(self_new_order)
    logger.debug(
        strm(
            "list searched:",
            output_shape[:-2] + [dot_dim_arg] + [output_shape[-1]],
        )
    )
    arg_new_order = tuple(
        arg.axn(j)
        for j in output_shape[:-2] + [dot_dim_arg] + [output_shape[-1]]
        if j in arg.dimlabels
    )
    logger.debug(
        strm("arg_new_order", arg_new_order, "arg.data.shape", arg.data.shape)
    )
    arg_formult_data = arg.data.transpose(arg_new_order)
    logger.debug(
        strm("new orders for self and arg are:", self_new_order, arg_new_order)
    )
    # {{{ insert singleton dims where needed
    for j, thisdim in enumerate(output_shape[:-2] + mult_dims[0]):
        if thisdim not in self.dimlabels and thisdim != rename_redundant[1]:
            self_formult_data = np.expand_dims(self_formult_data, j)
    for j, thisdim in enumerate(output_shape[:-2] + mult_dims[1]):
        if thisdim not in arg.dimlabels:
            arg_formult_data = np.expand_dims(arg_formult_data, j)
    # }}}
    logger.debug(
        strm(
            "raw nddata about to be multiplied are",
            self_formult_data.shape,
            arg_formult_data.shape,
        )
    )
    time_matmulstart = time.time()
    self.data = np.matmul(self_formult_data, arg_formult_data)
    # {{{ remove singleton dimensions that we added
    if "XXX_ADDED_XXX" in output_shape:
        this_slice = tuple(
            0 if j == "XXX_ADDED_XXX" else slice(None, None)
            for j in output_shape
        )
        self.data = self.data[tuple(this_slice)]
        output_shape = [j for j in output_shape if j != "XXX_ADDED_XXX"]
    # }}}
    time_matmulfinish = time.time()
    self.dimlabels = output_shape
    # {{{ use the dictionaries to reconstruct the metadata
    self.axis_coords = []
    self.axis_coords_units = []
    self.axis_coords_error = []
    for j, thisdim in enumerate(self.dimlabels):
        if thisdim == rename_redundant[1]:
            thisdim = rename_redundant[0]
        self.axis_coords.append(axis_coords_dict[thisdim])
        self.axis_coords_units.append(axis_units_dict[thisdim])
        self.axis_coords_error.append(axis_coords_error_dict[thisdim])
    # }}}
    logger.debug(strm("self.data.shape", self.data.shape))
    logger.debug(strm("ndshape", self.shape))
    time_dotfinish = time.time()
    logger.debug(
        f"initial took {(time_matmulstart-time_dotstart)*1e6:0.2f} μs, matmul {(time_matmulfinish-time_matmulstart)*1e6:0.2f} μs, finish {(time_dotfinish-time_matmulfinish)*1e6:0.2f} μs "
    )
    return self


# @profile
def matmul(self, arg):
    assert type(arg) is type(
        self
    ), "currently matrix multiplication only allowed if both are nddata"
    logger.debug("about to call dot")
    return self.C.dot(arg)
