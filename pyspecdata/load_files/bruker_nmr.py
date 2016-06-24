def series(filename, dimname=''):
    "For opening Bruker ser files"
    #{{{ Bruker 2D
    v = bruker_load_acqu(filename)
    v2 = bruker_load_acqu(filename,whichdim='2')
    td2 = int(v['TD'])
    rg = bruker_det_rg(float(v['RG']))
    td1 = int(v2['TD'])
    td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
    fp = open(filename+'ser','rb')
    data = fp.read()
    data = array(struct.unpack('>%di'%(len(data)/4),data),
            dtype='complex128')
    data = data[0::2]+1j*data[1::2]
    data /= rg
    mydimsizes = [td1,td2_zf/2]
    mydimnames = [dimname]+['t2']
    #print 'DEBUG: data going to nddata =',data
    try:
        data = nddata(data,mydimsizes,mydimnames)
    except:
        size_it_should_be = array(mydimsizes).prod()
        if size_it_should_be > len(data):
            zero_filled_data = zeros(size_it_should_be)
            zero_filled_data[0:len(data)] = data
            data = nddata(zero_filled_data,mydimsizes,mydimnames)
        else:
            new_guess = len(data)/(td2_zf/2)
            print lsafen("WARNING!, chopping the length of the data to fit the specified td1 of ",td1,"points!\n(specified ",zip(mydimnames,mydimsizes),' td2_zf=%d)'%td2_zf)
            #td2_zf_new = 2**ceil(log(td2)/log(2))
            #mydimsizes[1] = td2_zf_new
            #size_it_might_be = array(mydimsizes).prod()
            #print "maybe this works:",size_it_might_be == len(data)
            data = data[0:size_it_should_be]
            data = nddata(data,mydimsizes,mydimnames)
            #raise CustomError("found td1=",td1,"for",filename,"which I don't think is right, because the product of the dimensions",zip(mydimnames,mydimsizes),'=',size_it_should_be,'does not equal the length of the data',len(data),'I think that it should be',len(data)/(td2_zf/2))
    #print 'DEBUG: data straight from nddata =',data
    data = data['t2',0:td2/2] # now, chop out their zero filling
    t2axis = 1./v['SW_h']*r_[1:td2/2+1]
    t1axis = r_[0:td1]
    mylabels = [t1axis]+[t2axis]
    data.labels(mydimnames,mylabels)
    shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the first order phase shift
    data.circshift('t2',shiftpoints)
    data.set_units('t2','s')
    data.set_units('digital')
    data.other_info['title'] = bruker_load_title(filename)
    #print 'DEBUG 2: data from bruker file =',data
    #}}}
