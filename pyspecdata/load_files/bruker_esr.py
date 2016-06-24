def winepr(filename, dimname=''):
    "For opening WinEPR files"
    fp = open(filename+'.spc','rb')
    data = fp.read()
    data = array(
            struct.unpack('<%df'%(len(data)/4),data),
            dtype='double')
    v = winepr_load_acqu(filename)
    xpoints = v['RES']
    rg = v['RRG']
    data /= rg
    modulation = v['RMA']
    #data /= modulation
    try:
        data /= v['JNS'] # divide by number of scans
    except:
        pass
    #data /= v['MP'] # divide by power <-- weird, don't do this!
    ypoints = len(data)/xpoints
    if ypoints>1:
        if ypoints != v['REY']:
            raise CustomError('I thought REY was the indirect dim, guess not')
        if dimname=='':
            dimname = v['JEY']
        data = nddata(data,[ypoints,xpoints],[dimname,b0])
    else:
        data = nddata(data,[xpoints],[b0])
    xlabels = linspace(v['HCF']-v['HSW']/2.,v['HCF']+v['HSW']/2.,xpoints)
    if len(data.dimlabels)>1:
        yaxis = r_[0:v['REY']]
        if dimname == 'mw-power-sweep':
            yaxis *= v['MPS']
            yaxis += v['XYLB'] # the starting attenuation
            yaxis = 10**(-yaxis/10.) # convert to linear power
            yaxis *= v['MP']/yaxis[0] # the initial power
            yaxis *= 1e-3 # convert from mW to W
            data.rename('mw-power-sweep','power')
            dimname = 'power'
        data.labels([dimname,b0],[yaxis,xlabels])
        data.reorder([b0,dimname])
    else:
        data.labels([b0],[xlabels])
    return data
