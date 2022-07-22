from ..general_functions import *
if not inside_sphinx():
    from pylab import r_
def linear_shear(self,along_axis,propto_axis,shear_amnt,zero_fill = True):
    r"the linear shear -- see `self.shear` for documentation"
    if zero_fill:
        self.extend_for_shear(along_axis,propto_axis,shear_amnt)
    u_propto  = self.getaxis(propto_axis)# we will iterate along the propto dimension
    du_propto = check_ascending_axis(u_propto,additional_message="In order to implement a linearly interpolated skew with propto dimension %s"%propto_axis)
    du_along = check_ascending_axis(self.getaxis(along_axis),additional_message="In order to implement a linearly interpolated skew with along dimension %s"%along_axis)
    data_idx   = [s_[:]]*len(self.dimlabels)
    along_axn  = self.axn(along_axis)
    propto_axn = self.axn(propto_axis)
    # {{{ in the slice expression below, if propto_axn comes before along_axn, then the resulting slice will have reduced dimesionality
    if along_axn > propto_axn:
        along_axn -= 1
    # }}}
    # {{{ loop along the propto dimension (u_propto) and apply the shift
    for j in range(len(u_propto)):
        data_idx[propto_axn] = j
        integral_roll = u_propto[j] * shear_amnt // du_along # need to re-calculate
        #                                         for each slice, since it can vary
        decimal_remainder = u_propto[j] * shear_amnt / du_along - integral_roll
        integral_roll = int32(integral_roll)
        temp_floor = roll(self.data[data_idx], integral_roll,
                along_axn)
        temp_ceil = roll(self.data[data_idx], integral_roll + 1,
                along_axn)
        self.data[data_idx] = (1.-decimal_remainder) * temp_floor + decimal_remainder * temp_ceil
    # }}}
    return self
