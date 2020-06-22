from ..general_functions import *
if not inside_sphinx():
    from pylab import r_

def inhomog_coords(self, direct_dim, indirect_dim, tolerance = 1e-5,
        method = 'linear', plot_name = None, fl = None,
        debug_kwargs = {}):
    r'''Apply the "inhomogeneity transform," which rotates the data by :math:`45^{\circ}`, and then mirrors the portion with :math:`t_2<0` in order to transform from a :math:`(t_1,t_2)` coordinate system to a :math:`(t_{inh},t_{homog})` coordinate system.
    
    Parameters
    ----------
    direct_dim : str
        Label of the direct dimension (typically :math:`t_2`)
    indirect_dim : str
        Label of the indirect dimension (typically :math:`t_1`)
    method : 'linear', 'fourier'
        The interpolation method used to rotate the data and to mirror the data.
        **Note** currently, both use a fourier-based mirroring method.
    plot_name : str
        the base name for the plots that are generated
    fl : figlist_var
    debug_kwargs : dict
        with keys:

        :`correct_overlap`: if False, doesn't correct for the overlap error that occurs during mirroring
        '''
    correct_overlap = process_kwargs([('correct_overlap',True)],debug_kwargs)
    du = check_ascending_axis(self.getaxis(direct_dim),tolerance,"In order to perform the required rotations")
    print("first point along t2",self.getaxis(direct_dim)[0])
    assert abs(self.getaxis(direct_dim)[0]/du) < tolerance, "the "+direct_dim+"axis needs to start exactly at 0 -- you could try doing a sinc interpolation like this:\n\t" + "self.ft(direct_dim).ft_clear_startpoints(direct_dim,t=0,f='current').ift(direct_dim)"
    #{{{ and rotate!
    if fl is not None: fl.next(plot_name+', rotated')
    #{{{ calculate rotation parameters
    a = sqrt(2)-1. # = tan(pi/4./2.)
    b = -1./sqrt(2) # = -sin(pi/4.)
    #}}}
    if correct_overlap: self[direct_dim,0] *= 0.5 # correct the overlap artifact
    self.shear(direct_dim, indirect_dim, a, method = method)
    self.shear(indirect_dim, direct_dim, b, method = method)
    self.shear(direct_dim, indirect_dim, a, method = method)
    if fl is not None: fl.image(self)
    #}}}
    if fl is not None: fl.next(r'left side')
    if method == 'fourier':
        self.register_axis({direct_dim:0})
    left_data = self.copy()
    left_data[direct_dim,lambda x: x>0] = 0
    if fl is not None: fl.image(left_data)
    if fl is not None: fl.next(plot_name+': left side -- flipped')
    left_data.ft([indirect_dim,direct_dim])
    left_data.run(conj)
    left_data.ift([indirect_dim,direct_dim])
    if fl is not None: fl.image(left_data)
    if fl is not None: fl.next(plot_name+': right side')
    self[direct_dim,lambda x: x<0] = 0
    if fl is not None: fl.image(self)
    newdata = self[direct_dim,lambda x: x>=0]
    newdata += left_data[direct_dim,lambda x: x>=0]
    self.data = newdata.data
    self.setaxis(direct_dim,newdata.getaxis(direct_dim))
    print("At END: first point along t2",self.getaxis(direct_dim)[0])
    return self
