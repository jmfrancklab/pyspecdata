from .core import *
from scipy.interpolate import griddata
from matplotlib.tri import Triangulation, UniformTriRefiner
r'Helper functions for processing HFSS field data'

z_shift = 0.0

def gaussian_over(power_density,power_scale = None,
        f = 95e9,
        unit_scaling = 1e-3,# length units of mm
        w0 = 2.0,
        z0 = 148., # the z offset -- i.e. the z-value where the waist is
        ):
    r'''In order to compare the results of a simulation to a Gaussian beam,
    generate the power corresponding to a Gaussian beam over the same
    coordinates used in `power_density`.
    Following [1]_ (Chpt 2), the parameters used here are the *only* parameters needed
    to determine the fields for a Gaussian beam in free space.
    We use the condensed form [2]_:

    .. math::

        P=\frac{k^2}{2 a^2 \pi} \exp\left(\frac{-k^2 r^2}{2 a^2}\right)

    where

    .. math::

        a\equiv \sqrt{\left(\frac{1}{2} k w_0\right)^2+(z-z_0)^2}

    Parameters
    ==========
    power_density : nddata
        The data over which we want to generate the Gaussian beam profile.
        Only only the coordinate labels are used.
    power_scale : None or double
        If given, this sets the maximum power (at the center of the waist) of
        the resulting Gaussian beam.
        This is the same as replacing :math:`\frac{k^2}{2 a^2 \pi}` above with
        power_scale.
    f : double
        The frequency of the radiation in Hz.
    unit_scaling : double
        Determines the units that `w0` and `z0` are specified in.
        ``1e-3`` gives mm.
    z0 : double
        :math:`z_0`: The :math:`z` position of the beam waist (the
        narrowest part of the Gaussian beam).
    w0 : double
        :math:`w_0`: The beam waist radius, which is the radius of the
        beam at `z0`.
    
    References
    ----------
    .. [1] Goldsmith, P. Quasioptical Systems: Gaussian Beam Quasioptical Propogation and Applications; Wiley-IEEE Press, 1998.
    .. [2] Franck, J.M and Freed J.H. Resonator Design for High :math:`B_1` fields at 95 GHz *In Preparation* 2016
    '''
    c = 2.99792458e8
    k = 2*pi*f/c # in SI
    k *= unit_scaling
    a0 = 0.5*k*w0 # just a at z=0, which is just the raleigh length
    if 'z' in power_density.dimlabels and 'x' in power_density.dimlabels:
        a = power_density.fromaxis('z',lambda z: sqrt(a0**2+(z-z0)**2))
        r = power_density.fromaxis('x')
    elif power_density.get_prop('w_name') == 'z':
        z = power_density.get_prop('w_slice_position')/1e-3 # convert to mm
        print "pulled slice position as",z
        a = sqrt(a0**2+(z-z0)**2)
        r = power_density.fromaxis(['x','y'],lambda x,y:
                sqrt((x)**2+(y)**2))
    else:
        raise ValueError("I don't know how to interpret this data!!!")
    power = -k**2*r**2/2/a**2
    power.run(exp)
    power *= k**2/(2*a**2*pi)
    if power_scale is not None:
        power *= power_scale / (k**2/(2*a0**2*pi))
    return power
def w_index(u_index,v_index):
    'given indeces (0,1,2 for x,y,z) for u and v, return the index for the remaining axis -- *i.e.* "w"'
    return list(set([0,1,2])^set([u_index,v_index]))[0]
def det_w_name(u_label,v_label):
    'given labels x,y, and/or z for u and v, return the label of the remaining axis -- *i.e.* "w"'
    return ({'x','y','z'}^{u_label,v_label}).pop() # since there's only one element in the list, pop gives it
def construct_axes_from_positions(positions):
    r"""Take the array of `positions`, and use it to generate axes.
    Here, :math:`w` is the dimension along which the slice is taken, and :math:`u` and :math:`v` are the other two dimensions (typically in alphabetical order). 

    Parameters
    ----------
    positions : array
        an Nx3 array list of positions

    Returns
    -------
    w_axis_pos : double
        The position of the cross-sectional slice.
    u_index : int
        An index (0,1,2 for :math:`x`,:math:`y`,:math:`z`) that gives the
        identity of the first non-cross-sectional slice.
    u_axis : array
        A 1-D array giving the axis corresponding to `u_index`.
    v_index : int
        An index (0,1,2 for :math:`x`,:math:`y`,:math:`z`) that gives the
        identity of the second non-cross-sectional slice.
    v_axis : array
        A 1-D array giving the axis corresponding to `v_index`.
    """
    #{{{ construct the axes
    list_of_axes = [unique(positions[:,j]) for j in range(0,3)]
    list_of_indeces = [0,1,2]
    for j in range(0,len(list_of_axes)):
        if len(list_of_axes[j]) == 1:
            w_axis_pos = list_of_axes.pop(j)[0]
            list_of_indeces.pop(j)
            break
    [u_axis,v_axis] = list_of_axes
    [u_index,v_index] = list_of_indeces
    return w_axis_pos,u_index,v_index,u_axis,v_axis
    #}}}
def load_hfss_vectors(filename,show_valid = False,verbose = True):
    r"""Used by the main function, :func:`load_fields <pyspecdata.acert_hfss.load_fields>`, to open an HFSS .fld file that contains vector data:

    Returns
    =======
    w_axis_pos : double
        The position of the cross-sectional slice.
    u_index: uint
    u_axis: 
        :math:`u` and :math:`v` are the two axes (for, *e.g.*, an :math:`x-z` plane)
    v_index: uint
    v_axis:
        :math:`u` and :math:`v` are the two axes (for, *e.g.*, an :math:`x-z` plane)
    data: array
        a ``len(u_axis)`` :math:`\times` ``len(v_axis)`` :math:`\times 3` array
    """
    fp = open(filename,'r')
    data = fp.readlines()
    fp.close()
    header = data[0:2]
    if verbose: print "the second line tells me what's in the file -- x,y,z, followed by the 3 positions of vector data and the final part tells me what it is"
    if verbose: print 'header is',r'\begin{verbatim}', header, r'\end{verbatim}'
    data = data[2:]
    positions = empty((len(data),3),dtype = 'double')
    vec_vals = empty((len(data),3),dtype = 'double')
    for j,line in enumerate(data):
        vals = [double(x) for x in line.strip().split(" ") if len(x)>0]
        positions[j,:] = vals[0:3]
        vec_vals[j,:] = vals[3:]
    w_axis_pos,u_index,v_index,u_axis,v_axis = construct_axes_from_positions(positions)
    if verbose: print "the position along the $w$-axis (where $w$ is the cross-section dimension) is ",w_axis_pos
    #{{{ show the datapoints and which are valid
    if show_valid:
        fl.next('show valid values')
        thismask = isnan(vec_vals[:,2])
        fl.plot(positions[:,0][thismask],positions[:,v_index][thismask],'r.')
        fl.plot(positions[:,0][~thismask],positions[:,v_index][~thismask],'b.')
        xlabel('x')
        ylabel('z')
        axes().set_aspect('equal', 'datalim')
    #}}}
    data = empty((len(u_axis),len(v_axis),3),dtype = 'double')
    for j in range(0,positions.shape[0]):
        u_i = where(u_axis == positions[j,u_index])[0][0]
        v_i = where(v_axis == positions[j,v_index])[0][0]
        data[u_i,v_i,:] = vec_vals[j,:]
    return w_axis_pos,u_index,u_axis,v_index,v_axis,data
def load_hfss_scalar(filename,verbose = False):
    r"""Used by the main function, :func:`load_fields <pyspecdata.acert_hfss.load_fields>`.  Like load_hfss_vectors, but returns a scalar, so the third (spatial) dimension is dropped"""
    fp = open(filename,'r')
    data = fp.readlines()
    fp.close()
    header = data[0:2]
    if verbose: print "the second line tells me what's in the file -- x,y,z, followed by the 3 positions of vector data and the final part tells me what it is"
    if verbose: print 'header is',r'\begin{verbatim}',header,r'\end{verbatim}'
    data = data[2:]
    positions = empty((len(data),3),dtype = 'double')
    value = empty(len(data),dtype = 'double')
    for j,line in enumerate(data):
        vals = [double(x) for x in line.strip().split(" ") if len(x)>0]
        positions[j,:] = vals[0:3]
        value[j] = vals[3]
    w_axis_pos,u_index,v_index,u_axis,v_axis = construct_axes_from_positions(positions)
    data = empty((len(u_axis),len(v_axis)),dtype = 'double')
    for j in range(0,positions.shape[0]):
        u_i = where(u_axis == positions[j,u_index])[0][0]
        v_i = where(v_axis == positions[j,v_index])[0][0]
        data[u_i,v_i] = value[j]
    return w_axis_pos,u_index,u_axis,v_index,v_axis,data
def load_fields_numpyversion(basename = 'mesh_only_140129',slicename = 'xz',fieldname = 'poynting',has_mesh = True):
    r"""this is obsolete -- was the original version of load_fields, where I was not using nddata

    Returns
    =======
    u_index: int
    v_index: int
    U: array 
        array of the :math:`u`-coordinate labels -- matches the dimensions of `data`
    V: array
        array of the :math:`v`-coordinate labels -- matches the dimensions of `data`
    data: array
        array giving the 2D data stored in the .fld file
    model_data: array
        array giving the 2D model data stored in the .fld file
    """
    def load_hfss_complex_vectors(filename):
        w_axis_pos,u_index,u_axis,v_index,v_axis,data_real = load_hfss_vectors(filename%"real")
        w_axis_pos,u_index,u_axis,v_index,v_axis,data_imag = load_hfss_vectors(filename%"imag")
        data = data_real + 1j*data_imag
        return u_index,u_axis,v_index,v_axis,data

    u_index,u_axis,v_index,v_axis,data = load_hfss_complex_vectors(getDATADIR()+'HFSS_jmf356/'+fieldname+'_%s_'+slicename+'_'+basename+'.fld') # format string where "real" and "imag" go
    _,_,_,_,_,temp = load_hfss_scalar(getDATADIR()+'HFSS_jmf356/real_dielectric_'+slicename+'_'+basename+'.fld')
    temp[~isfinite(temp)] = -1
    model_data = zeros_like(temp,dtype = [('mesh','d'),('dielectric','d')])
    model_data['dielectric'] = temp
    if has_mesh:
        _,_,_,_,_,temp = load_hfss_scalar(getDATADIR()+'HFSS_jmf356/mesh_domain_'+slicename+'_'+basename+'.fld')
        temp[~isfinite(temp)] = -1
        model_data['mesh'] = temp
    U = u_axis[:,None]*ones(len(v_axis))[None,:]
    V = v_axis[None,:]*ones(len(u_axis))[:,None]
    return u_index,v_index,U,V,data,model_data
def load_fields(basename = 'mesh_only_140129',slicename = 'xz',fieldname = 'poynting',has_mesh = True):
    r"""Finds the fld files named by `basename`, `slicename`, and `fieldname`, and loads them appropriately.

    It assumes that there are the following files:

        #. ``[fieldname]_real_[slicename]_[basename].fld`` where ``[fieldname]`` is something like "poynting" or "E," while "slicename" is something like "xy" or "xy_zoom"
        #. same format, with ``real`` :math:`\rightarrow` ``imag``, for the imaginary part
        #. ``[fieldname]_real`` :math:`\rightarrow` ``real_dielectric`` for the real dielectric (used to determine the model) 
        #. (if `has_mesh` is `True`) ``[fieldname]_real`` :math:`\rightarrow` ``mesh_domain`` for the mesh domain (in general, this can be any conductor domain that's called 'mesh domain' and is also used to determine the model, since the dielectric isn't usually sufficient)

    It does this with relatively few lines of code, relying on
    :func:`load_hfss_vectors <pyspecdata.acert_hfss.load_hfss_vectors>`
    and
    :func:`load_hfss_scalar <pyspecdata.acert_hfss.load_hfss_scalar>`
    to do the heavy lifting,
    but it should be the function typically used.

    Returns
    =======
    data: nddata
        array giving the 2D data stored in the .fld file, with the following properties set:
            :w_name: The name of the slice dimension (*e.g.*, if `data` has two dimensions `x` and `y`, then this is :math:`z`.
            :w_slice_position: The position of the slice along the dimension given by `w_name`
        Aside from the two spacial dimensions, it has a third dimension, called
        `vector`, which gives the three vector components (:math:`x`, :math:`y`,
        and :math:`z`) at the particular coordinate.
        Sets the `name` property to ``basename+' '+slicename``.
    model_data: nddata
        In the same format as data, except that it has no `vector` dimension.
        Rather, the data is a structured numpy array with the field names
        `mesh` and `dielectric`, which give the location of the mesh and
        dielectric model elements.
        Sets the `name` property to ``basename+' '+slicename``.
    """
    def cartesian_bynumber(a):
        return ['x','y','z'][a]
    def load_hfss_complex_vectors(filename):
        w_axis_pos,u_index,u_axis,v_index,v_axis,data_real = load_hfss_vectors(filename%"real")
        w_axis_pos,u_index,u_axis,v_index,v_axis,data_imag = load_hfss_vectors(filename%"imag")
        data = data_real + 1j*data_imag
        return w_axis_pos,u_index,u_axis,v_index,v_axis,data

    w_axis_pos,u_index,u_axis,v_index,v_axis,data = load_hfss_complex_vectors(getDATADIR()+'HFSS_jmf356/'+fieldname+'_%s_'+slicename+'_'+basename+'.fld') # format string where "real" and "imag" go
    _,_,_,_,_,temp = load_hfss_scalar(getDATADIR()+'HFSS_jmf356/real_dielectric_'+slicename+'_'+basename+'.fld')
    temp[~isfinite(temp)] = -1
    model_data = zeros_like(temp,dtype = [('mesh','d'),('dielectric','d')])
    model_data['dielectric'] = temp
    if has_mesh:
        _,_,_,_,_,temp = load_hfss_scalar(getDATADIR()+'HFSS_jmf356/mesh_domain_'+slicename+'_'+basename+'.fld')
        temp[~isfinite(temp)] = -1
        model_data['mesh'] = temp
    uv_names = map(cartesian_bynumber,[u_index,v_index])
    return nddata(data,data.shape,uv_names+['vector']).set_prop(dict(
        w_name = cartesian_bynumber(w_index(u_index,v_index)),
        w_slice_position = w_axis_pos)
        ).labels(uv_names,[u_axis,v_axis]).name(basename+' '+slicename), nddata(
                model_data,model_data.shape,uv_names).set_prop(dict(
                    w_name = cartesian_bynumber(w_index(u_index,v_index)),
                    w_slice_position = w_axis_pos)
                    ).labels(uv_names,[u_axis.copy(),v_axis.copy()]).name(basename+' '+slicename)
def contour_power(S_data, model_data,
        rotate=False, direction_of_power=None,
        number_of_arrows=1, figsize=None, amount_above=0.1,
        power_scale=None, fl=None, z_shift=None, residual=True,
        z_limit=None,
        **kwargs):
    r"""Plot the magnitude of the power flow, determined from the Poynting vector (:math:`\vec{S}`), along with a contour intended to indicate a Gaussian beam radius.
    Plots either the HFSS field or a `residual` with a Gaussian beam (the residual is the default -- see below).

    It does this by calling :func:`load_fields <pyspecdata.acert_hfss.load_fields>`, placing the result into an appropriate nddata function, and then plotting the power density using an image function, and a contour using an nddata :func:`contour <pyspecdata.nddata.contour>` method.

    Parameters
    ----------
    S_data : nddata
        The Poynting vector data, in a format as generated by load_fields.
    model_data : nddata
        The model data, in a format as generated by load_fields.
    residual : bool
        Show the residual with a Gaussian beam determined by :func:`gaussian_over <pyspecdata.acert_hfss.gaussian_over>`
    z_shift : double
        (in mm) Set the :math:`z` value `z_shift` as the origin on the plot
    rotate : boolean
        rotate the plot by 90 degrees relative to its default
        orientation
    direction_of_power : {0,1,2,'x','y','z'}
        by default, it determines power flow along :math:`u \times v`, but this
        allows you to override by setting the direction of power to x, y, z, or
        the corresponding integer.
    power_scale : None or double
        optionally used to scale the overall power
    fl : figlist_var
        The figure list to plot to.  Use the `S_data.name()` as the figure name.
    kwargs : dict
        the remaining keyword arguments are passed on to
        :func:`gaussian_over <pyspecdata.acert_hfss.gaussian_over>`

    Returns
    -------
    Nothing, unless the power is flowing parallel to the plane, in which case,
    it returns the max power.
    """
    if fl is None:
        raise ValueError('please just pass a figure list (fl)')
    for thisaxis in S_data.dimlabels:
        if thisaxis != 'vector':
            S_data.set_units(thisaxis,'mm')
            model_data.set_units(thisaxis,'mm')
            S_data.setaxis(thisaxis,lambda x: x/1e-3)
            model_data.setaxis(thisaxis,lambda x: x/1e-3)
    if z_shift is not None:
        if 'z' in S_data.dimlabels:
            S_data.setaxis('z',lambda x: x-z_shift)
            model_data.setaxis('z',lambda x: x-z_shift)
        else:
            raise ValueError("You're trying to set z_shift, but the axes here are"+' and '.join(self.dimlabels))
    if rotate:
        S_data.reorder(S_data.dimlabels[1])
        model_data.reorder(model_data.dimlabels[1])
        if figsize is None:
            fl.next(S_data.name(),figsize = (5,6))
    else:
        if figsize is None:
            fl.next(S_data.name(),figsize = (5,12))
    if figsize is not None:
            fl.next(S_data.name(),figsize = figsize)
    #{{{ calculate the normalized power density
    normalize_byslice = False # i.e. assume power is flowing
    #                   through the plane, so we want to
    #                   normalize against the *overall* max
    if direction_of_power is None:
        direction_of_power_str = det_w_name(*tuple(model_data.dimlabels))# use model data, since S_data also include the vector dimension
    else:
        if type(direction_of_power) is str:
            direction_of_power_str = direction_of_power
        else:
            direction_of_power_str = ['x','y','z'][direction_of_power]
        if direction_of_power_str != det_w_name(*tuple(model_data.dimlabels)):
            normalize_byslice = True
    direction_of_power = ['x','y','z'].index(direction_of_power_str)
    power_density = S_data['vector',direction_of_power].runcopy(real) #  taking
    #                            just the perpendicular component
    #}}}
    max_power = power_density.copy()
    max_power[lambda x: ~isfinite(x)] = 0
    for j in power_density.dimlabels:
        if j != direction_of_power_str:
            max_power.run(max,j)
    normalized_power_density = power_density / max_power
    overall_max_power = None
    if normalize_byslice: # power flows parallel to the plot
        overall_max_power = max_power.data.max()
        if power_scale is None:
            power_scale = overall_max_power
    if power_scale is not None:
        overall_max_power = power_scale
    if residual:
        power_density_residual = power_density - gaussian_over(power_density,power_scale = overall_max_power,**kwargs)
        power_density_residual *= 100. / overall_max_power
        forplot = power_density_residual
    else:
        forplot = power_density
    if z_limit is not None and 'z' in forplot.dimlabels:
        print "applying z limit",z_limit
        forplot = forplot['z':(z_limit,)]
    fl.image(forplot,x_first = True)
    if z_limit is not None and 'z' in forplot.dimlabels:
        normalized_power_density = normalized_power_density['z':(z_limit,)]
    normalized_power_density.contour(colors='w', alpha=0.75,
            labels=False, levels=[exp(-2.0)], linewidths=3)
    if overall_max_power is not None:
        textstr = 'max power {:g}'.format(overall_max_power)
        if residual:
            textstr += '\nvalues give % residual'
        text(0.3*amount_above, 0.0 + 0.3*amount_above,# these say from the left and from the bottom of the figure -- 1.05 allows space for the title
                textstr,
                horizontalalignment = 'left',
                verticalalignment = 'bottom',
                size = 'x-small',
                transform = gca().transAxes)
    if len(max_power.data.shape) < 1:
        text(1.0, 1.0 - amount_above,# these say from the left and from the bottom of the figure -- 1.05 allows space for the title
                'max power along %s (normalized) is %g'%(['x','y','z'][direction_of_power],max_power.data),
                horizontalalignment = 'right',
                verticalalignment = 'bottom',
                size = 'x-small',
                transform = gca().transAxes)
    if rotate:
        gca().invert_yaxis()
        fig = gcf()
        fig.subplots_adjust(left = 0.22, right = 0.98, top = 0.825, bottom =0.17)# to allow space for labels
        axes().set_aspect('equal', 'datalim')
    else:
        axes().set_aspect('equal', 'datalim')
        fig = gcf()
        fig.subplots_adjust(left = 0.3, right = 0.95, top = 0.825, bottom =0.075)# to allow space for labels
    if normalize_byslice:
        return overall_max_power
    else:
        return
