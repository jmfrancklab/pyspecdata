#from .core import *
#from scipy.interpolate import griddata
#from matplotlib.tri import Triangulation, UniformTriRefiner
r'Helper functions for processing HFSS field data'

z_shift = 0.0

def gaussian_over(power_density,power_scale = None):
    r'''In order to compare the results of a simulation to a Gaussian beam,
    generate the power corresponding to a Gaussian beam over the same
    coordinates used in `power_density`.

    Parameters
    ==========
    power_density : nddata
        The data over which we want to generate the Gaussian beam profile.
        Only only the coordinate labels are used.
    power_scale : None or double
        If given, this sets the maximum power (at the center of the waist) of
        the resulting Gaussian beam.
    
    '''
    k = 2*pi/3e-3 # in SI
    k *= 1e-3 # length units of mm
    w0 = 2.0
    z0 = 148. # the z offset -- i.e. the z-value where the waist is
    a0 = 0.5*k*w0 # just a at z=0, which is just the raleigh length
    if 'z' in power_density.dimlabels and 'x' in power_density.dimlabels:
        a = power_density.fromaxis('z',lambda z: sqrt(a0**2+(z-z0)**2))
        r = power_density.fromaxis('x')
    elif 'x' in power_density.dimlabels and 'y' in power_density.dimlabels:
        a = 1.
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
def construct_axes_from_positions(positions):
    'take the array of `positions`, and use it to generate 1D axes for u and v'
    #{{{ construct the axes
    list_of_axes = [unique(positions[:,j]) for j in range(0,3)]
    list_of_indeces = [0,1,2]
    for j in range(0,len(list_of_axes)):
        if len(list_of_axes[j]) == 1:
            list_of_axes.pop(j)
            list_of_indeces.pop(j)
            break
    [u_axis,v_axis] = list_of_axes
    [u_index,v_index] = list_of_indeces
    return u_index,v_index,u_axis,v_axis
    #}}}
def load_hfss_vectors(filename,show_valid = False,verbose = False):
    r"""open an HFSS .fld file that contains vector data:

    Returns
    =======
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
    u_index,v_index,u_axis,v_axis = construct_axes_from_positions(positions)
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
    return u_index,u_axis,v_index,v_axis,data
def load_hfss_scalar(filename,verbose = False):
    ("like load_hfss_vectors, but returns a scalar, so the third"
    "(spatial) dimension is dropped")
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
    u_index,v_index,u_axis,v_axis = construct_axes_from_positions(positions)
    data = empty((len(u_axis),len(v_axis)),dtype = 'double')
    for j in range(0,positions.shape[0]):
        u_i = where(u_axis == positions[j,u_index])[0][0]
        v_i = where(v_axis == positions[j,v_index])[0][0]
        data[u_i,v_i] = value[j]
    return u_index,u_axis,v_index,v_axis,data
def load_fields(basename = 'mesh_only_140129',slicename = 'xz',fieldname = 'poynting',has_mesh = True):
    r"""Finds the fld files named by `basename`, `slicename`, and `fieldname`, and loads them appropriately.

    It assumes that there is a file called ``[fieldname]_real_[slicename]_[basename].fld``, as well as one for the imaginary part, as well as one for the real dielectric (used to determine the model) and (if `has_mesh` is `True`), one for the mesh domain (in general, this can be any conductor domain that's called 'mesh domain' and is also used to determine the model, since the dielectric isn't usually sufficient)

    It does this with relatively few lines of code, relying on
    :func:`load_hfss_vectors <pyspecdata.acert_hfss.load_hfss_vectors>`
    and
    :func:`load_hfss_scalar <pyspecdata.acert_hfss.load_hfss_scalar>`
    to do the heavy lifting.

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
        u_index,u_axis,v_index,v_axis,data_real = load_hfss_vectors(filename%"real")
        u_index,u_axis,v_index,v_axis,data_imag = load_hfss_vectors(filename%"imag")
        data = data_real + 1j*data_imag
        return u_index,u_axis,v_index,v_axis,data

    u_index,u_axis,v_index,v_axis,data = load_hfss_complex_vectors(getDATADIR()+'HFSS_jmf356/'+fieldname+'_%s_'+slicename+'_'+basename+'.fld') # format string where "real" and "imag" go
    junk1,junk2,junk3,junk4,temp = load_hfss_scalar(getDATADIR()+'HFSS_jmf356/real_dielectric_'+slicename+'_'+basename+'.fld')
    temp[~isfinite(temp)] = -1
    model_data = zeros_like(temp,dtype = [('mesh','d'),('dielectric','d')])
    model_data['dielectric'] = temp
    if has_mesh:
        junk1,junk2,junk3,junk4,temp = load_hfss_scalar(getDATADIR()+'HFSS_jmf356/mesh_domain_'+slicename+'_'+basename+'.fld')
        temp[~isfinite(temp)] = -1
        model_data['mesh'] = temp
    U = u_axis[:,None]*ones(len(v_axis))[None,:]
    V = v_axis[None,:]*ones(len(u_axis))[:,None]
    return u_index,v_index,U,V,data,model_data
def contour_power(basename=None, slicename=None,
        rotate=False, direction_of_power=None,
        number_of_arrows=1, figsize=None, amount_above=0.1,
        power_scale = None,
        fl = None,
        **kwargs):
    r"""Plot the magnitude of the power flow, determined from the Poynting vector (:math:`\vec{S}`), along with a contour intended to indicate a Gaussian beam radius.

    It does this by calling :func:`load_fields <pyspecdata.acert_hfss.load_fields>`, placing the result into an appropriate nddata function, and then plotting the power density using an image function, and a contour using an nddata :func:`contour <pyspecdata.nddata.contour>` method.

    Parameters
    ----------
    z_shift : double
        Set the :math:`z` value `z_shift` as the origin on the plot
    rotate : boolean
        rotate the plot by 90 degrees relative to its default
        orientation
    direction_of_power : integer
        by default, it determines power flow along
        :math:`u \times v`,
        but this allows you to override
    power_scale : None or double
        optionally used to scale the overall power

    Returns
    -------
    Nothing, unless the power is flowing parallel to the plane, in which case,
    it returns the max power.
    """
    if fl is None:
        raise ValueError('please just pass a figure list (fl)')
    u_index,v_index,U,V,S_data,model_data = load_fields(basename = basename,slicename = slicename,fieldname = 'poynting',**kwargs)
    U /= 1e-3
    V /= 1e-3
    V -= z_shift
    if rotate:
        u_index,v_index = v_index,u_index
        U,V = V.T,U.T
        S_data = S_data.transpose([1,0,2])
        model_data = model_data.T
        if figsize is None:
            fl.next(basename+' '+slicename,figsize = (5,6))
    else:
        if figsize is None:
            fl.next(basename+' '+slicename,figsize = (5,12))
    if figsize is not None:
        fl.next(basename+' '+slicename,figsize = figsize)
    #{{{ calculate the normalized power density
    normalize_byslice = False # i.e. assume power is flowing
    #                   through the plane, so we want to
    #                   normalize against the *overall* max
    if direction_of_power is None:
        direction_of_power = w_index(u_index,v_index)
    else:
        if direction_of_power != w_index(u_index,v_index):
            normalize_byslice = True
    power_density = real(S_data[:,:,direction_of_power])#  taking
    #                            just the perpendicular component
    #}}}
    axis_labels = map(lambda x: ['x','y','z'][x],[u_index,v_index])
    power_density = nddata(power_density, shape(power_density),
            axis_labels).labels(axis_labels,[U[:,0],V[0,:]])
    max_power = power_density.copy()
    max_power[lambda x: ~isfinite(x)] = 0
    direction_of_power_str = ['x','y','z'][direction_of_power]
    for j in axis_labels:
        if j != direction_of_power_str:
            max_power.run(max,j)
    normalized_power_density = power_density / max_power
    if normalize_byslice: # power flows parallel to the plot
        overall_max_power = max_power.data.max()
        if power_scale is None:
            power_scale = overall_max_power
    if power_scale is not None:
        overall_max_power = power_scale
    power_density_residual = power_density - gaussian_over(power_density,power_scale = overall_max_power)
    power_density_residual *= 100. / overall_max_power
    fl.image(power_density_residual,x_first = True)
    normalized_power_density.contour(colors='w', alpha=0.75,
            labels=False, levels=[exp(-2.0)], linewidths=3)
    if len(max_power.data.shape) < 1:
        text(1.0, 1.0 + amount_above,# these say from the left and from the bottom of the figure -- 1.05 allows space for the title
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
