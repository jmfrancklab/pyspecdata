<<<<<<< HEAD
import os
#import termios
import socket
import time
from struct import *
from pylab import *
from socket import socket
import scipy.io as mio
=======
from matlablike import *
import os
import re
import socket
import struct
import time

>>>>>>> public
DEFAULTCOM = 4

#{{{ the underlying gpib class, which has older instruments integrated (separate them later
class gpib():
<<<<<<< HEAD
    def __init__(self,ip = '192.168.0.2',timeout = 0.1):
        # Switch for OS X
        self.flags = {}
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
        self.socket.settimeout(timeout)
        self.socket.connect((ip, 1234))
        
        # Set mode as CONTROLLER
        self.socket.write('++mode 1'+"\r")
        # what does this do?
        self.socket.write('++ifc'+"\r")
        # Turn off read-after-write to avoid "Query Unterminated" errors
        self.socket.write('++auto 0'+"\r")
        # Read timeout is 500 msec
        sock.send("++read_tmo_ms 500\n")
        # ethernet gpib code (arb_eth.py) does the following, but I'm not going to, because it's not like before
        ## Do not append CR or LF to GPIB data
        #sock.send("++eos 3\n")
        ## Assert EOI with last byte to indicate end of data
        #sock.send("++eoi 1\n")
        self.socket.write('++eoi 0'+"\r")
        self.caddr = -1
        self.socket.write("++ver\r")
        versionstring = self.socket.readline()
        if versionstring[0:8]=='Prologix':
            print 'connected to: ',versionstring
        else:
            print 'Error! can\'t find prologix on COM%d'%(usbnumber+1)
            raise
    def close(self):
        self.serial.write('++clr'+"\r")
        self.serial.write('++loc'+"\r")
        self.serial.close()
        
    def setaddr(self,addr):
        if(self.caddr != addr):
            self.serial.write('++addr '+str(addr)+"\r")
            time.sleep(0.1)
            self.caddr=addr
    def readline(self,addr):    
        self.setaddr(addr)
        self.serial.write('++read 10'+"\r")
        return self.serial.readline()[:-2] # get rid of dos newline
    def read(self,addr,numberofbytes = None):
        self.setaddr(addr)
        self.serial.write('++read eoi'+"\r")
        if numberofbytes:
            return self.serial.read(numberofbytes)
        else:
            return self.serial.readline()
    def write(self,addr,gpibstr):
        self.setaddr(addr)
        self.serial.write(gpibstr+"\r")
    def respond(self,addr,gpibstr,printstr = '%s'):
        self.serial.write('++auto 1'+"\r") # because sometimes, it's just talking at me nonstop
        self.serial.flush() # so go ahead and flush all that garbage
        self.write(addr,gpibstr) # now write whatever we want
        retval = printstr % self.read(addr) # and instantly read the id string
        self.serial.write('++auto 0'+"\r") # now set back into the normal mode, where I need to request a response
        return retval
        
    #{{{ Functions for Newer Tek Scope
    def tek_query_var(self,addr,varname):
        self.write(addr,varname+'?')
        temp = self.read(addr)
        temp = temp[len(varname)+2:-1] # remove initial space and trailing \n
        if temp[0]=='\"':
            return temp[1:-1]
        else:
            return double(temp)
    def tek_get_curve(self,addr):
        y_unit = self.tek_query_var(addr,'WFMP:YUN')
        y_mult = self.tek_query_var(addr,'WFMP:YMU')
        y_offset = self.tek_query_var(addr,'WFMP:YOF')
        dx = self.tek_query_var(addr,'WFMP:XIN')
        x_unit = self.tek_query_var(addr,'WFMP:XUN')
        #print y_mult,y_unit,y_offset,dx,x_unit
        self.write(addr,'CURV?')
        self.serial.write('++addr '+str(addr)+"\r")
        time.sleep(0.1)
        self.serial.write('++read eoi'+"\r")
        header_string = self.serial.read(8)
        print "'"+header_string+"'\n"
        print "reading length of length: "+header_string[-1]
        curve_length = self.serial.read(int(header_string[-1]))
        print "reading curve of length: "+curve_length
        x = header_string[0:int(curve_length)]*dx
        return (x_unit,
                y_unit,
                x,
                y_offset+y_mult*array(
                unpack(
                    '%sb'%curve_length,
                    self.serial.read(int(curve_length))
                    )))
    #}}}
    
    #{{{ Functions for HP 54110D Digitizing Oscilloscope
    def hp_get_curve(self,addr):
        # Acquire waveform
        self.write(addr,":acquire:type normal")
        self.write(addr,":digitize CHANNEL2")
        self.write(addr,":waveform:source 2")
        self.write(addr,":waveform:format ascii")
        
        self.write(addr,":waveform:data?");
        
        self.write(addr,"++read eoi");
        wfrm = [int(x) for x in self.serial.readlines()]
        
        # Acquire preamble
        self.write(addr,":waveform:preamble?")
        self.write(addr,"++read eoi");
        pra = self.serial.readline().split(",")
        print 'pra=\'',pra,'\''
        try:
            format = int(pra[0])
            type = int(pra[1])
            points = int(pra[2])
            count = int(pra[3])
        
            xinc = float(pra[4])
            xorig = float(pra[5])
            xref = int(pra[6])
        
            yinc = float(pra[7])
            yorig = float(pra[8])
            yref = int(pra[9])
        
        except IndexError:
            print "Bad preamble recieved"
            exit(1)
        
        if points != len(wfrm):
            print "WARNING: Received less points than specified in the preamble"
        
        x = ((r_[0:len(wfrm)]-xref)*xinc)+xorig
        y = ((array(wfrm)-yref)*yinc)+yorig
        
        # FIXME: No idea what x_unit and y_unit are for. They just get stowed
        #        in the matlab file so for now it's okay. /eazg
        return (1,1,x,y)
    #}}}
#}}}

#{{{ a copy of the above, for an ethernet-based controller
class gpib_eth():
    def __init__(self,address="192.168.0.100",port=1234):
        # Switch for OS X
        self.flags = {}
        self.port = socket()
        self.port.connect((address,port))
        self.port.send('++mode 1'+"\r")
        self.port.send('++ifc'+"\r")
        self.port.send('++auto 0'+"\r")
        self.port.send('++eoi 0'+"\r")
        self.port.send("++ver\r")
        versionstring = self.port.recv(1000)
        self.caddr = -1
        if versionstring[0:8]=='Prologix':
            print 'connected to: ',versionstring
        else:
            print 'Error! can\'t find prologix on %s:%d'%(address,port)
            raise
    def close(self):
        self.port.close()
        
    def setaddr(self,addr):
        if(self.caddr != addr):
            self.port.send('++addr '+str(addr)+"\r")
            self.caddr=addr
    def readandchop(self,addr): # unique to the ethernet one
        retval = self.port.recv(1024) # get rid of dos newline
        while (retval[-1] == '\r') or (retval[-1] == '\n'): # there should be a function for this (i.e. chop, etc)!
            retval = retval[:-1]
        return retval
    def readline(self,addr):    
        self.setaddr(addr)
        self.port.send('++read 10'+"\r")
        return self.readandchop(addr)
    def read(self,addr):
        self.setaddr(addr)
        self.port.send('++read eoi'+"\r")
        return self.readandchop(addr)
    def write(self,addr,gpibstr):
        self.setaddr(addr)
        self.port.send(gpibstr+"\r")
    def respond(self,addr,gpibstr,printstr):
        self.write(addr,gpibstr)
        print printstr % self.read(addr)
        
    #{{{ Functions for Newer Tek Scope
    def tek_query_var(self,addr,varname):
=======
    def __init__(self,ip = '127.0.0.1',timeout = 4.0,buffer_len = 512,sock = 1234,version = True):
        # Switch for OS X
        self.buffer_len = buffer_len
        self.flags = {}
        self._timeout = timeout
        #self.socket.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
        self._sock = sock
        self._ip = ip
        self.caddr = -1 # current gpib address
        self._connect()
        # ethernet gpib code (arb_eth.py) does the following, but I'm not going to, because it's not like before
        #self.socket.send("++eot_enable 0\n") # Do not append user-specified character to the end
        if version == True:
            self.ver()
        return
    def _connect(self):
        r'This both reopens the socket, and also resets the gpib address to self.caddr'
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
        self.socket.settimeout(self._timeout)
        try:
            self.socket.connect((self._ip, self._sock))
        except:
            extratext = ''
            if self._ip == '127.0.0.1':
                extratext = '\nmake sure that the ssh tunnels are up and running!!!'
            raise RuntimeError('connection to the prologix connector at ip address %s:%d failed!'%(self._ip,self._sock)+extratext)
        self.socket_connected = True
        self.setaddr(self.caddr,force = True)
        self.socket.send('++mode 1'+"\n") # Set mode as CONTROLLER
        self.socket.send('++ifc'+"\n") # this does a low-level GPIB thing to grab control of the GPIB bus
        self.socket.send('++auto 0'+"\n") # this buffers the data (until ++read eoi) in the prologix, rather than automatically spitting it out
        self.set_timeout(1000e-3) # Read timeout is 500 msec
        self.socket.send("++eos 3\n") # Do not append CR or LF to GPIB data
        self.socket.send("++eoi 1\n") # Assert EOI with last byte to indicate end of data -- this allows me to buffer the data in the prologix and read out with ++read eoi
        return
    def disconnect_socket(self):
        self.socket.close()
        del self.socket
        self.socket_connected = False
        return
    def ver(self):
        if not self.socket_connected: self._connect()
        self.socket.send("++ver\n")
        versionstring = self.socket.recv(self.buffer_len)
        if versionstring[0:8]=='Prologix':
            print 'connected to: ',versionstring,'on IP',self._ip,':',self._sock
        else:
            raise RuntimeError('Error! can\'t find prologix on %s:%d'%(self._ip,self._sock))
        return
    def close(self):
        raise ValueError('I changed this --> you no longer call instance.close() and instead, you say del instance')
    def __del__(self):
        r'close the connection to the Prologix adapter'
        if not self.socket_connected: self._connect()
        self.socket.send('++clr'+"\n")
        self.socket.send('++loc'+"\n")
        self.socket.send('++rst'+"\n") # Power-on reset of the controller
        self.disconnect_socket()
        return
    def setaddr(self,addr,force = False):
        r'set the current GPIB address to "addr"'
        if not self.socket_connected: self._connect()
        if (force and addr > 0) or (self.caddr != addr):
            self.socket.send('++addr '+str(addr)+"\n")
            self.caddr=addr
        return
    def set_timeout(self,timeout):
        if timeout is None:
            timeout = 500e-3
        timeout = int(1e3*timeout)
        if timeout>3000:
            print 'timeout of greater than 3000ms not allowed'
            timeout = 3000
        if not self.socket_connected: self._connect()
        self.socket.send("++read_tmo_ms %d\n"%timeout)
        print 'set prologix timeout to %d ms'%timeout
        return
    def readline(self,addr = None):
        r'read a line of text coming from the current gpib address'
        if addr is not None:
            self.setaddr(addr)
        if not self.socket_connected: self._connect()
        self.socket.send('++read 10'+"\n")
        retval = self.socket.recv(self.buffer_len)
        return retval
    def serial_poll(self,addr = None):
        if addr is not None:
            self.setaddr(addr)
        retval = self.respond('++spoll'+"\n")
        return (unpackbits(uint8(retval))[::-1] == True).tolist()
        #return map(lambda x: bool(x&retval),
        #        map(lambda x: 1<<x,
        #            range(0,8)))
    def oldread(self,addr = None,timeout = None):
        r'read from the current gpib address until the end of input'
        if not self.socket_connected: self._connect()
        if timeout is not None:
            oldtimeout = self.socket.gettimeout()
            self.socket.settimeout(timeout)
        if addr is not None:
            self.setaddr(addr)
        self.socket.send('++read eoi'+"\n")
        #time.sleep(0.1)
        mystr = self.socket.recv(self.buffer_len)
        #try:
        #except:
        #    raise ValueError('I tried to read from the GPIB device, but it timed out')
        retstr = mystr
        #while len(mystr) == self.buffer_len:
        while True:
            try:
                mystr = self.socket.recv(self.buffer_len)
            except:
                break
            retstr = retstr + mystr
        if timeout is not None:
            self.socket.settimeout(oldtimeout)
        return retstr
    def read(self,addr = None,timeout = None,until = 'eoi',length = None,keep_trying = False,prologix_timeout = None,go_until_term = True,verbose = False): #Remember to change keep_trying back to False
        r'''read from the current gpib address
keyword args:
	until --> character (or eoi) to read up until
	timeout --> if set, overrides current socket timeout
	keep_trying --> keep trying to read until the socket times out
	prologix_timeout --> the timeout when the prologix talks to the instrument
	go_until_term --> if the terminating character does not appear at the end of the buffer, keep reading from the buffer until it does
'''
        if prologix_timeout is not None:
            self.set_timeout(prologix_timeout)
        until_is_char = False
        if len(until) == 1:
            until_chr = until
            until = '%d'%ord(until)
            until_is_char = True
        if length is not None:
            until_chr = ''
        if timeout is not None:
            if verbose: print 'resetting timeout for read to',timeout
            oldtimeout = self.socket.gettimeout()
            self.socket.settimeout(timeout)
        if addr is not None:
            self.setaddr(addr)
        if not self.socket_connected: self._connect()
        self.socket.send("++read %s \n"%until)
        #time.sleep(0.1)
        #try:
        mystr = self.socket.recv(self.buffer_len)
        #except:
        #    raise ValueError('I tried to read from the GPIB device, but it timed out')
        if verbose: print 'len of data',len(mystr),'of buffer',self.buffer_len
        retstr = mystr
        while len(mystr) == self.buffer_len and length is None:
                try:
                    if verbose: print('I\'m about to read the stuff from the buffer')
                    mystr = self.socket.recv(self.buffer_len)
                    if verbose: print('mystr: ' + mystr)
                    if verbose: print 'reloaded buffer, len =',len(mystr)
                except:
                    break
                retstr = retstr + mystr
        if (until_is_char and go_until_term) or (length is not None and len(retstr) >= self.buffer_len):
            while retstr[-1] != until_chr:
                try:
                    mystr = self.socket.recv(self.buffer_len)
                    if verbose: print 'read %d more seeking end, total length is currently %d'%(len(mystr),len(retstr))
                except:
                    raise ValueError('Failed while attempting to read more until I hit the term character,'+repr(until_chr))
                retstr = retstr + mystr
                if length is not None and len(retstr) == length:
                    break
        if keep_trying:
            print 'keep trying to grab data, even after it seems to be gone'
            while True:
                try:
                    mystr = self.socket.recv(self.buffer_len)
                    if verbose: print 'reloaded buffer, len =',len(mystr)
                except:
                    break
                retstr = retstr + mystr
        if prologix_timeout is not None:
            self.set_timeout(500e-3)
        if verbose: print 'read() data:',retstr
        if verbose: print 'len of data:',len(retstr)
        return retstr
    def write(self,gpibstr,addr = None,verbose = False):
        r'write gpibstr to the current gpib address'
        if addr is not None:
            self.setaddr(addr)
        string = gpibstr + '\n'
        if verbose:
            print 'Sending:' + string
        if not self.socket_connected: self._connect()
        bytes_out = self.socket.send(gpibstr+"\n")
        if bytes_out != len(string):
            raise ValueError('Network is busy...')
        return
    def respond(self,gpibstr,printstr = '%s',addr = None,wait = False,timeout = None,until = '\n',length = None,prologix_timeout = None,verbose = False):
        'write gpibstr to the current gpib address, then print the response, optionally formatted according to the format given in printstr\nmore recently, I added functionality where it will wait (in 5 second increments if wait is set to true, or longer if you set it to a number)\nby default, it reads the response up until the newline (\\n)\nhowever, you can also set until = "eoi" for binary data, where applicable'
        if not self.socket_connected: self._connect()
        self.write(gpibstr,addr = None,verbose = verbose) # now write whatever we want
        if not wait:
            retval = printstr % self.read(addr = addr,timeout = timeout,until = until,length = length,prologix_timeout = prologix_timeout,verbose = verbose) # and instantly read the id string
        else:
            if wait is True:
                wait = 5
            try:
                retval = printstr % self.read(addr = addr,timeout = timeout,until = until,prologix_timeout = prologix_timeout,verbose = verbose) # and instantly read the id string
            except:
                print 'waiting for a response I know will come'
                time.sleep(wait)
                try:
                    retval = printstr % self.read(addr = addr,timeout = timeout,until = until,length = length,prologix_timeout = prologix_timeout,verbose = verbose) # and instantly read the id string
                except:
                    print 'still waiting for a response I know will come'
                    time.sleep(wait)
                    try:
                        retval = printstr % self.read(addr = addr,timeout = timeout,until = until,prologix_timeout = prologix_timeout,verbose = verbose) # and instantly read the id string
                    except:
                        print 'still waiting for a response I know will come, but now I am losing hope'
                        time.sleep(wait)
                        try:
                            retval = printstr % self.read(addr = addr,timeout = timeout,until = until,prologix_timeout = prologix_timeout,verbose = verbose) # and instantly read the id string
                        except:
                            raise ValueError('I waited a very long time for the GPIB to respond, but I lost hope')
        retval = retval.strip()
        try:
            new_retval = int(retval)
        except:
            try:
                new_retval = double(retval)
            except:
                try:
                    retval_list = retval.split(',')
                except:
                    return retval #if it's not a list, return
                if retval_list[-1] == '':
                    retval_list.pop(-1)
                try:
                    new_retval = array(map(int,retval_list))
                except:
                    try:
                        new_retval = array(map(double,retval_list))
                    except:
                        new_retval = retval
        return new_retval
        
    #{{{ Functions for Newer Tek Scope
    def tek_query_var(self,addr,varname):
        r'an old function for talkint to Devin\'s Tek oscilloscope'
>>>>>>> public
        self.write(addr,varname+'?')
        temp = self.read(addr)
        temp = temp[len(varname)+2:-1] # remove initial space and trailing \n
        if temp[0]=='\"':
            return temp[1:-1]
        else:
            return double(temp)
    def tek_get_curve(self,addr):
<<<<<<< HEAD
=======
        r'an old function for talkint to Devin\'s Tek oscilloscope'
>>>>>>> public
        y_unit = self.tek_query_var(addr,'WFMP:YUN')
        y_mult = self.tek_query_var(addr,'WFMP:YMU')
        y_offset = self.tek_query_var(addr,'WFMP:YOF')
        dx = self.tek_query_var(addr,'WFMP:XIN')
        x_unit = self.tek_query_var(addr,'WFMP:XUN')
        #print y_mult,y_unit,y_offset,dx,x_unit
        self.write(addr,'CURV?')
<<<<<<< HEAD
        self.serial.write('++addr '+str(addr)+"\r")
        time.sleep(0.1)
        self.serial.write('++read eoi'+"\r")
=======
        self.serial.write('++addr '+str(addr)+"\n")
        time.sleep(0.1)
        self.serial.write('++read eoi'+"\n")
>>>>>>> public
        header_string = self.serial.read(8)
        print "'"+header_string+"'\n"
        print "reading length of length: "+header_string[-1]
        curve_length = self.serial.read(int(header_string[-1]))
        print "reading curve of length: "+curve_length
        x = header_string[0:int(curve_length)]*dx
        return (x_unit,
                y_unit,
                x,
                y_offset+y_mult*array(
<<<<<<< HEAD
                unpack(
=======
                struct.unpack(
>>>>>>> public
                    '%sb'%curve_length,
                    self.serial.read(int(curve_length))
                    )))
    #}}}
    
    #{{{ Functions for HP 54110D Digitizing Oscilloscope
    def hp_get_curve(self,addr):
        # Acquire waveform
        self.write(addr,":acquire:type normal")
        self.write(addr,":digitize CHANNEL2")
        self.write(addr,":waveform:source 2")
        self.write(addr,":waveform:format ascii")
        
        self.write(addr,":waveform:data?");
        
        self.write(addr,"++read eoi");
        wfrm = [int(x) for x in self.serial.readlines()]
        
        # Acquire preamble
        self.write(addr,":waveform:preamble?")
        self.write(addr,"++read eoi");
        pra = self.serial.readline().split(",")
        print 'pra=\'',pra,'\''
        try:
            format = int(pra[0])
            type = int(pra[1])
            points = int(pra[2])
            count = int(pra[3])
        
            xinc = float(pra[4])
            xorig = float(pra[5])
            xref = int(pra[6])
        
            yinc = float(pra[7])
            yorig = float(pra[8])
            yref = int(pra[9])
        
        except IndexError:
            print "Bad preamble recieved"
            exit(1)
        
        if points != len(wfrm):
            print "WARNING: Received less points than specified in the preamble"
        
        x = ((r_[0:len(wfrm)]-xref)*xinc)+xorig
        y = ((array(wfrm)-yref)*yinc)+yorig
        
        # FIXME: No idea what x_unit and y_unit are for. They just get stowed
        #        in the matlab file so for now it's okay. /eazg
        return (1,1,x,y)
    #}}}
#}}}

#{{{ this section gives wrappers for specific instruments
class eip_powermeter ():
<<<<<<< HEAD
    def __init__(self,comport = DEFAULTCOM,gpibaddress=19):
        self.g = gpib(comport - 1) # on port 4 with the current connector
        self.gpibaddress=gpibaddress
        self.g.write(self.gpibaddress,'R5')# output at a lower resolution for faster sampling
        self.g.write(self.gpibaddress,'PR')# undocumented, outputs just power reading
        #self.g.write(self.gpibaddress,'DP')# turn off the display
        self.g.write(self.gpibaddress,'PA')# make sure the power meter is on
        #self.g.write(self.gpibaddress,'FA')# set to "fast mode"??
        self.g.write(self.gpibaddress,'HP')# set the hold off
        self.g.write(self.gpibaddress,'RA')# stay in "data output mode"
    def read_power(self):
        self.g.write(self.gpibaddress,'RS')# "reset" which apparently takes a reading
        try:
            retval = float(self.g.readline(self.gpibaddress))
        except:
            retval = -999.9
        counter = 0
        while (counter < 4) & (retval == -999.9):
            #print 'reading...'
            #self.g.write(self.gpibaddress,'RS')# "reset" which apparently takes a reading
            tempstr = self.g.readline(self.gpibaddress)
            if len(tempstr)>0:
                retval = float(tempstr)
            else:
                retval = -999.9
            counter += 1
            print '/',
            time.sleep(1e-4)
        if retval == -999.9:
            print 'failed to read a power!'
        return retval
    def close(self):
=======
    def __init__(self,gpibaddress=19,sock = 1234,ip = '127.0.0.1'):
        self.g = gpib(ip = ip,sock = sock,version = False)
        self.gpibaddress = gpibaddress
        self.g.setaddr(self.gpibaddress)
        self.g.write('R5') # set resolution to 5 figures
        self.g.write('PA') # set to power and frequency mode
        self.g.write('PR') # set to power reading mode
    def read_power(self):
        power = self.g.readline() # read oout power
        power.strip() # clean up string
        power = float(power) # convert to float
        return power
    def __del__(self):
>>>>>>> public
        self.g.write(self.gpibaddress,'DA')# turn on the display
        self.g.write(self.gpibaddress,'HP')# if we don't do this, the display freezes
        self.g.write(self.gpibaddress,'RP')# no longer output only mode
        self.g.write(self.gpibaddress,'FP')# turn off "fast mode"??
        self.g.write(self.gpibaddress,'R0')# switch back to high res
<<<<<<< HEAD
        self.g.close()
=======
        del self.g
        
>>>>>>> public
class hp8672a():
    def __init__(self,comport = DEFAULTCOM,gpibaddress = 3):
        self.g = gpib(comport - 1) # on port 4 with the current connector
        self.gpibaddress = gpibaddress
        self.stepsize = 0.5e6 # this is a lie, but it's used as a default by wobbandmin
    def set_frequency(self,frequency):
        self.g.write(self.gpibaddress,'P%08dZ0'%int(round(frequency*1e-3)))# just use the 10 GHz setting, and fill out all the other decimal places with zeros
        return
    def close(self):
<<<<<<< HEAD
        self.g.close()
#}}}

#{{{ here define wrapper function that automatically identify different types of instruments
def powermeter(comport = DEFAULTCOM,gpibaddressrange = range(1,21)):
    for gpibaddress in gpibaddressrange:
        print "trying address",gpibaddress
        g = gpib(comport - 1,timeout = 0.1) # use a timeout of half a second, otherwise this is painfully slow
        idstring = g.respond(gpibaddress,'ID')
        g.close()
        if idstring[0:4] == 'GIGA':
            return gigatronics_powermeter(comport = comport,gpibaddress = gpibaddress)
        elif idstring.find('E0') > 0: #because I can't get EIP to return an id string right now
            return eip_powermeter(comport = comport,gpibaddress = gpibaddress)
#}}}

#{{{ this section gives classes for specific instruments

class gigatronics_powermeter ():

    def __init__(self,comport = DEFAULTCOM,gpibaddress=15):
        self.g = gpib(comport - 1) # on port 4 with the current connector
        self.gpibaddress=gpibaddress

        idstring = self.g.respond(self.gpibaddress,'ID') # Check ID command
        if idstring[0:4] == 'GIGA':
            print 'idstring is',idstring
            self.g.write(self.gpibaddress,'TR3')        # Set Free Run Trigger Mode
            self.g.write(self.gpibaddress,'LG')         # Set Log units in dBm
=======
        raise ValueError('I changed this --> you no longer call instance.close() and instead, you say del instance')
    def __del__(self):
        del self.g
#}}}

#{{{ here define wrapper function that automatically identify different types of instruments
def powermeter(comport = DEFAULTCOM,gpibaddressrange = range(1,36)):
    ip = '127.0.0.1'
    sock = 1234
    for gpibaddress in gpibaddressrange:
        print "trying address",gpibaddress
        try:
            g = gpib(ip = ip,sock = sock,timeout = .5) # use a timeout of half a second, otherwise this is painfully slow
            #idstring = g.respond(gpibaddress,'ID')
            idstring = g.respond('ID')
            print('I got the id string:', idstring)
            if idstring[0:4] == 'GIGA':
                print('\nfound GIGA at %i!!!!!!!!!!!!!!!\n' %gpibaddress)
                return gigatronics_powermeter(comport = comport,gpibaddress = gpibaddress)
            elif idstring.find('E0') > 0: #because I can't get EIP to return an id string right now
                print('not ok')
                #return eip_powermeter(comport = comport,gpibaddress = gpibaddress)
            
            del g
        except:
            print('Couldn\'t find anything here')
#}}}

#{{{ this section gives classes for specific instruments
class agilent():
    #def __init__(self,sock = 5025,timeout = 0.5,buffer_len = 4096,ip = '127.0.0.1'):
    def __init__(self,sock = 5025,timeout = 10.0,buffer_len = 4096,ip = '127.0.0.1'):
        self.buffer_len = buffer_len
        self.flags = {}
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM, socket.IPPROTO_TCP)
        self.socket.settimeout(timeout)
        self.socket.connect((ip, 5025))
        self.socket.send('*IDN?\n')
        versionstring = self.socket.recv(self.buffer_len)
        print 'Connected to',versionstring
    def command_send(self,cmd):
        self.socket.send(cmd)
    def command(self,cmd):
        self.socket.send(cmd)
        return self.socket.recv(self.buffer_len)
    def readvoltage(self):
        self.socket.send(':MEAsure:VMIN? CHANnel4\n')
        voltage = self.socket.recv(self.buffer_len)
        return voltage
    def setvoltage(self,ch1V,ch2V = None):
        if ch2V == None:
            ch2V = ch1V
        self.command_send(':CHANnel1:SCALe %.2fV\n'%ch1V)
        self.command_send(':CHANnel2:SCALe %.2fV\n'%ch2V)
        return
    def clear(self):
        self.socket.send(':DISPlay:CLEar\n')
    def position(self, position):
        self.socket.send(':TIMebase:POSition %0.3e\n'%position)
    def read(self,terminator = None,length = None,verbose = False):
        if (terminator is not None) and (length is not None):
            raise ValueError('Both terminator and length can\'t be set to a value!!!')
        try:
            mystr = self.socket.recv(self.buffer_len)
        except:
            raise ValueError('Scope timed out when trying to read')
        
        retstr = mystr
        while len(mystr) == self.buffer_len:
            if verbose: print('len(mystr) == self.buffer_len')
            try:
                mystr = self.socket.recv(self.buffer_len)
            except:
                break
            retstr = retstr + mystr
        if terminator is not None:
            while retstr[-1] != terminator:
                try:
                    mystr = self.socket.recv(self.buffer_len)
                except:
                    break
                retstr = retstr + mystr
        if length is not None:
            while len(retstr) < length:
                if verbose: print('len(retstr) < length')
                mystr = self.socket.recv(self.buffer_len)
                retstr = retstr + mystr
        return retstr
    def timebase(self, resolution):
        # sets timebase in s/div
        resolution *= 10
        self.socket.send(':TIMebase:RANGe %0.3e\n'%resolution)

    def waveform_complex(self, avg_pts = 1):
        data = self.waveform(1, avg_pts)
        data.data = data.data + 1j*self.waveform(2, avg_pts).data
        return data

    def waveform(self, channel, avg_pts = 1, preamble = None, get_preamble = False):
        avg_pts *= 1024
        save_timeout = self.socket.gettimeout()
        long_timeout = 1.4 * avg_pts / 1024
        self.socket.send(':WAVeform:FORMat WORD\n')
        self.socket.send(':WAVeform:BYTeorder LSBFirst\n')
        ydata = [0, 0, 0, 0, 0, 0, 0]
        self.socket.send(':MTESt:AVERage:COUNt ' + str(int(avg_pts)) + '\n')
        
        self.socket.send(':WAVeform:SOURce CHANnel' + str(channel) + '\n')
        if preamble == None:
            self.socket.send(':WAVeform:PREAmble?\n')
            if save_timeout < long_timeout:
                self.socket.settimeout(long_timeout)
            header = self.read(terminator = '\n')
            self.socket.settimeout(save_timeout)
        else:
            header = preamble
        
        mylist = header.split(',')
        points = int(mylist[2])
        xscale = float(mylist[4])
        xorigin = float(mylist[5])
        yscale = float(mylist[7])
        yorigin = float(mylist[8])
        resol = int(mylist[9])
        
        self.socket.send(':WAVeform:DATA?\n')
        if save_timeout < long_timeout:
           self.socket.settimeout(long_timeout)
        trace = self.read(length = 10+points*2)
        self.socket.settimeout(save_timeout)
        
        data = trace[10:points*2+10]
        ydata = fromstring(data,dtype='uint16') * yscale - resol * yscale + yorigin
        xdata = arange(0,xscale*points,xscale)[0:points]
        
        if get_preamble:
            return nddata(ydata,[-1],['t'],axis_coords = [xdata]), header
        else:
            return nddata(ydata,[-1],['t'],axis_coords = [xdata])
    def acquire(self,counts = 10000):
        '''setup sampling scope to acquire data with the digitize function'''
        #self.socket.send('*CLS\n') #clear status
        #self.socket.send(':DISPlay:CLEar\n') #clear the display and reset all associated measurements

        #self.socket.send(':TIMebase:RANGe 5e-7\n') #time base to 200 ns/div (range = 10*time/div) 
        #self.socket.send(':TIMebase:DElay 132.5e-9\n') #delay to zero
        #self.socket.send(':TIMebase:REFerence CENTer\n') #display ref. at center

        #self.socket.send(':CHANel1:PROBe 10') #probe attenuation to 10:1
        #self.socket.send(':CHANel1:RANGe 8\n') #vertical range 1.6 v full scale #gives undefined header
        #self.socket.send(':CHANel1:OFFSet 1.06250\n') #offset to -0.4
        #self.socket.send(':CHANel1:COUPling dc\n') #coupling to dc

        #self.socket.send(':TRIGger:SWEep NORmal\n') #normal triggering
        #self.socket.send(':TRIGger:LEVel -.626\n') #trigger level to -0.4
        #self.socket.send(':TRIGger:SLOPe POSitive\n') #trigger on positive slope
        
        #self.socket.send(':acquire:type normal\n') #normal acquisition
        #self.socket.send(':ACQuire:TYPE NORMal\n') 
        #self.socket.send(':ACQuire:TYPE AVERage\n') 
        #self.socket.send(':ACQuire:COMPlete 100\n')
        self.socket.send(':ACQuire:COUNt %i\n'%counts)
        self.xdata = self.get_preamble(1)

        self.socket.send(':WAVeform:FORMat WORD\n')
        self.socket.send(':WAVeform:BYTeorder LSBFirst\n')
        self.socket.send(':WAVeform:UNSigned 1\n') #Turn unsigned return off (0) or on (1)

    def digitize(self):
        '''acquire data on all channels and save in buffer'''
        self.socket.send(':DIGitize\n')
    def run(self):
        '''Equivalent to pressing run on agilent'''
        self.socket.send(':RUN\n')
    def Waveform(self,channel,do_get_preamble = False):
        '''Tim's waveform function modified from thomas\'
        this pulls the data from the agilent scope on a given channel'''
        self.socket.send(':WAVeform:SOURce CHANnel%i\n' %(channel))
        if do_get_preamble:
            self.xdata = self.get_preamble(channel)

        self.socket.send('WAVeform:DATA?\n')
        trace = self.read(length = 10+2*self.points)
        self.ydata = fromstring(trace[10:self.points*2+10],dtype='uint16') * self.yscale - self.resol * self.yscale + self.yorigin
        return self.xdata,self.ydata
    def get_preamble(self,channel):
        self.socket.send(':WAVeform:SOURce CHANnel%i\n' %(channel))
        self.socket.send(':WAVeform:PREamble?\n')
        header = self.read(terminator = '\n')

        mylist = header.split(',')
        self.points = int(mylist[2])
        self.xscale = float(mylist[4])
        self.xorigin = float(mylist[5])
        self.yscale = float(mylist[7])
        self.yorigin = float(mylist[8])
        self.resol = int(mylist[9])

        self.xdata = arange(0,self.xscale*self.points,self.xscale)[0:self.points]
        return self.xdata
    def Waveform_auto(self):
        self.digitize()
        try:
            self.xdata
        except:
            self.xdata = self.get_preamble
        xdata, real_y = self.Waveform(1)
        xdata, imag_y = self.Waveform(2)
        return nddata(real_y + 1j*imag_y,[-1],['t'],axis_coords = [xdata])
        
        
#{{{ instek scope
class instek ():
    def __init__(self,gpibaddress=5,sock = 1234,ip = '127.0.0.1'):
        self.g = gpib(ip = ip)
        self.gpibaddress = gpibaddress
        self.g.setaddr(self.gpibaddress)
        self.g.write('*idn?')
        idstring = self.g.read()
        print 'connected to ',idstring,'\n'
    
    def avg(self,channel):
        self.g.socket.send(':measure:source1 ch' + str(channel) + '\n')
        self.g.socket.send(':measure:average?\n')
        average = float(self.g.read())
        return average
    
    def waveform(self,channel):
        self.g.socket.send(':acquire' + str(channel) + ':memory?\n')
        trace = self.g.read()
        eoh = trace.find('#550000')
        header = trace[0:eoh+7]
        data = trace[eoh+7:]
        trace = [header, data[0:50000]]
        
        mylist = header.split(';')
        for i in range(len(mylist)):
            mylist[i] = mylist[i].split(',')
        yscale = float(mylist[6][1])
        xscale = float(mylist[10][1])
        xdata = arange(0,2*xscale,xscale/12500)
        ydata = fromstring(trace[1],dtype='>i2') / 25. * yscale
        print 'integral:',
        print sum(ydata)
        
        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata, ydata)
        
        ax.set_xlim(0, 2*xscale)
        ax.set_ylim(-4*yscale, 4*yscale)
        ax.set_xlabel('time [s]')
        ax.set_ylabel('voltage [v]')
        
        show()

        return trace

    def close(self):
        raise valueerror('i changed this --> you no longer call instance.close() and instead, you say del instance')
    def __del__(self):
        del self.g
#}}}
#{{{ gigatronics power meter
class gigatronics_powermeter ():
    def __init__(self,gpibaddress=15,sock = 1234,ip = '127.0.0.1'):
        self.g = gpib(ip = ip,sock = sock)
        self.gpibaddress = gpibaddress
        self.g.setaddr(self.gpibaddress)
        try:
            idstring = self.g.respond('ID')
        except:
            raise ValueError('The gigatronics powermeter not responding to the ID command --> maybe it\'s not on, not hooked up, or you\'re talking to the wrong GPIB address')

        if idstring[0:4] == 'GIGA':
            print 'idstring is',idstring
            self.g.write('TR3')        # Set Free Run Trigger Mode
            self.g.write('LG')         # Set Log units in dBm
>>>>>>> public
            #self.g.write(self.gpibaddress,'DD')         # Display Disable
        else:
            raise ValueError('Not a Gigatronics power meter, returned ID string %s'%idstring)
        
    def read_power(self):
        try:
            retval = float(self.g.readline(self.gpibaddress))
        except:
            retval = -999.9
        counter = 0
        while (counter < 4) & (retval == -999.9):
            #print 'reading...'
            #self.g.write(self.gpibaddress,'RS')# "reset" which apparently takes a reading
            tempstr = self.g.readline(self.gpibaddress)
            if len(tempstr)>0:
                retval = float(tempstr)
            else:
                retval = -999.9
            counter += 1
            print '/',
            time.sleep(1e-4)
        if retval == -999.9:
            print 'failed to read a power!'
        return retval
    
    def close(self):
<<<<<<< HEAD
        #self.g.write(self.gpibaddress,'DE')         # Display Enable
     
        
##        self.g.write(self.gpibaddress,'HP')# if we don't do this, the display freezes
##        self.g.write(self.gpibaddress,'RP')# no longer output only mode
##        self.g.write(self.gpibaddress,'FP')# turn off "fast mode"??
##        self.g.write(self.gpibaddress,'R0')# switch back to high res
        self.g.close()

=======
        raise ValueError('I changed this --> you no longer call instance.close() and instead, you say del instance')
    def __del__(self):
        del self.g
#}}}
#{{{ Stanford lock-in
class lockin_amp ():
    def __init__(self,gpibaddress=3,sock = 1234,ip = '127.0.0.1'):
        self.g = gpib(ip = ip,sock = sock)
        self.gpibaddress = gpibaddress
        self.g.setaddr(self.gpibaddress)
        try:
            idstring = self.g.respond('*IDN?')
        except:
            raise ValueError('The lockin amp is not responding to the ID command --> maybe it\'s not on, not hooked up, or you\'re talking to the wrong GPIB address')
        print 'Connected to ',idstring,'\n'
        self.front_display(1,'X')
        self.front_display(2,'Y')
        #{{{ set the filter very narrow, wait, and perform the auto offset
        self.wait_until_ready()
        self.filter(24)
        #{{{ perform auto sensitivity
        #self.auto_gain()
        self.wait_until_ready()
        #}}}
        ##{{{ auto offset, and redo sensitivity
        #self.auto_offset()
        #self.g.write('AGAN')
        #self.wait_until_ready()
        #print 'READY'
        ##}}}
    def serial_poll(self):
        return dict(zip(['SCN','IFC','ERR','LIA','MAV','ESB','SRQ','unused'],
            self.g.serial_poll()))
    def auto_gain(self):
        self.g.write('AGAN')
        self.wait_until_ready()
        return
    def expand(self,val = None,on = ['X','Y']):
        r'query the expand values on X and Y by default, and if val is given (can be a list), set them'
        retval = {}
        if val is not None:
            if type(val) is not list:
                val = [val]*len(on)
            val = [round(log(x)/log(10)) for x in val]
        for j,k in enumerate([1,2]):
            this_item = '%d,%d'%(k,self._convert_code(on[j],k))
            if val is not None:
                self.g.write('DEXP%s,%0.0f'%(this_item,val[j]))
            retval.update({
                on[j]:10**(int(self.g.respond('DEXP?%s'%this_item)))})
        return retval
    def auto_offset(self,on = ['X','Y']):
        r'perform auto offset (on X and Y by default)'
        for j,k in enumerate(on):
            self.g.write('AOFF %d,%d'%self._convert_code(on[j]))
            self.wait_until_ready()
        #}}}
    def wait_until_ready(self):
        while not self.serial_poll()['IFC']:
            print 'not ready yet, IFC is False'
            time.sleep(1e-3)
    def _convert_code(self,value,channel = None,units = 'Volts'):
        if value == 'R':
            value = 'R [%s]'%units
        codes_for_1 = {'X':0,
                'R [Volts]':1,
                'R [dBm]':2,
                'X noise':3,
                'AUX IN 1':4}
        codes_for_2 = {'Y':0,
                'theta':1,
                'Y noise [Volts]':2,
                'Y noise [dBm]':3,
                'AUX IN 2':4}
        need_to_know_channel = False
        if channel is None:
            need_to_know_channel = True
            if value in codes_for_1.keys():
                channel = 1
            elif value in codes_for_2.keys():
                channel = 2
            else:
                raise ValueError('You\'ve given me just the string value',value,'without specifying which channel you\'re talking about, which would be OK, except that neither channel 1',codes_for_1.keys(),'nor channel 2',codes_for_2.keys(),"can take this value")
        if channel == 1:
            codes = codes_for_1
        elif channel == 2:
            codes = codes_for_2
        if type(value) is str:
            try:
                thiscode = codes[value]
            except:
                raise ValueError('The value %s'%value+'needs to be one of'+','.join(codes.keys())+'when this is called for channel %d'%channel)
            if need_to_know_channel:
                return channel,thiscode
            else:
                return thiscode
        else:
            if need_to_know_channel:
                raise ValueError("You can't pass a list without passing the channel number")
            return [k for k,v in codes.iteritems() if v == value]
    def front_display(self,channel,value = None,units = 'Volts'):
        r'if value is not given, just query the value of the front display for the given channel\nif value is given, then set the front display to that value'
        if value is not None:
            self.g.write('DDEF%d,%d'%(channel,
                self._convert_code(value,channel,units = units)
                ))
        return self._convert_code(int(self.g.respond('DDEF?%d'%(channel))),channel)
    def close(self):
        raise ValueError('I changed this --> you no longer call instance.close() and instead, you say del instance')
    def __del__(self):
        del self.g
    def set_amp(self,sens):
        r'Change Sensitivity'
        self.g.write('SENS %0.0f' % (sens))
    def set_phase(self,phas):
        r'Set Phase to phas (must be between -179.99 and 180.00)'
        self.g.write('PHAS %0.2f' % (phas))
    def auto_phase(self):
        r'Same as pressing Shift-Phase. Adjusts reference phase so that the current measurement has a Y value of zero and an X value equal to the signal magnitude, R'
        self.g.write('APHS')
    def _conv_dboct_to_internal(self,value,back = False):
        r'convert to and from db/oct and the internal representation'
        if type(value) is str or type(value) is int:
            value = double(value)
        if back:
            return value*6
        else:
            retval = round(value/6)
            if retval > 4:
                retval = 4
            return retval
    def filter(self,*arg):
        r'set or query the time constant filter slope, either call as filter(k) to set the slope k in dB/oct'
        if len(arg) == 1:
            self.g.write('OFSL%d'%self._conv_dboct_to_internal(arg[0]))
        self.wait_until_ready()
        query_filter = self.g.respond('OFSL?',wait = True)
        return self._conv_dboct_to_internal(query_filter,back = True)
    def amp_out(self,voltage = None,aux = 1):
        r'''query or set the output amplitude on the auxilliary output,
        either call as amp_out(V) to set to voltage V
        or as amp_out() to query the current voltage'''
        if voltage is not None:
            self.g.write('AUXO%d,%f0.3'%(aux,voltage))
        return double(self.g.respond('AUXO?%d'%(aux)))
    def set_harmonic(self,mode = None):
        if mode is not None:
            mode -= 1
            self.g.write('HARM%d'%mode)
        return int(self.g.respond('HARM?'))+1
    def set_wavetype(self,wave):
        self.g.write('')
        return
    def read(self):
        return double(self.g.respond('OUTP?1')) + 1j * double(self.g.respond('OUTP?2'))
    def out_data(self):
        return self.g.respond('OUTR?1').strip()
    def _set_time_const(self,n = None):
        r'''Query time constant or set time constant
        n value     sample rate     n value     sample rate
        0           100us           9           3s
        1           300us           10          10s
        2           1ms             11          30s
        3           3ms             12          100s
        4           10ms            13          300s
        5           30ms            14          1ks
        6           100ms           15          3ks
        7           300ms           16          10ks
        8           1s              17          30ks'''
                   
        if n == None:
            return self.g.respond('OFLT?')
        else:
            n = str(int(n))
            self.g.write('OFLT ' + n)
    def _set_sample_rate(self,f = None):
        r'''Query sample rate or set the sample rate.
        set_sample_rate(f) where f is the sample rate in Hz.
        n value     sample rate
        0           62.5 mHz (16 sec)
        1           125 mHz (8 sec)
        2           250 mHz (4 sec)
        3           500 mHz (2 sec)
        4           1 Hz
        5           2 Hz
        6           4 Hz
        7           8 Hz
        8           16 Hz
        9           32 Hz
        10          64 Hz
        11          128 Hz
        12          256 Hz
        13          512 Hz
        14          Trigger'''
        if f == None:
            print 'Query sample rate '
            return self.g.respond('SRAT?')
        else:
            self.n = int(round(log2(float(f)) + 4))
            if self.n < 0 or self.n > 14:
                raise ValueError('Invalid input for f. int(round(log2(self.f) + 4)) must be an integer between 0 and 14')
            else:
                self.g.write('SRAT ' + self.n.__str__())
    def _set_gain(self,i = None):
        r'''Query gain or set gain
        i value     Sensitivity         i value     Sensitivity
        0           100nVrms/-127dBm    8           1mVrms/-47dBm
        1           300nVrms/-117dBm    9           3mVrms/-37dBm
        2           1uVrms/-107dBm      10          10mVrms/-27dBm
        3           3uVrms/-97dBm       11          30mVrms/-17dBm
        4           10uVrms/-87dBm      12          100mVrms/-7dBm
        5           30uVrms/-77dBm      13          300mVrms/+3dBm
        6           100uVrms/-67dBm     14          1Vrms/+13dBm
        7           300uVrms/-57dBm'''
        if i == None:
            return self.g.respond('SENS?')
        else:
            i = str(int(i))
            self.g.write('SENS ' + i)
    def _pause_storage(self):
        self.g.write('PAUS')
    def _reset_storage(self):
        self.g.write('REST')
    def retrieve_stored(self,verbose = False):
        print 'Start time: \n'
        if not hasattr(self,'time_start'):
            raise RuntimeError('Before you can retrieve stored, you must start_storage!!')
        print self.time_start
        # do not send base-level socket commands
        # instead, use the "self.g.write" and "self.g.respond" functions only
        self.g.write('PAUS') #Stop data acquisition
        number_of_points = int(self.g.respond('SPTS?')) #Query Number of Points
        if verbose: print('Number of points: ' + str(number_of_points))
        frequency = int(round(2**(float(self.g.respond('SRAT?'))-4))) #Query frequency and convert to Hz
        if verbose: print('Frequency:'+str(frequency))
        data_ch1 = self.g.respond('TRCA?' + '1,0,' + str(number_of_points),verbose = verbose) #Return Ch. 1
        data_ch2 = self.g.respond('TRCA?' + '2,0,' + str(number_of_points),verbose = verbose) #Return Ch. 2
        if verbose: print 'checkpoint 1'
        self.time_axis = r_[0:double(number_of_points)]/frequency+self.time_start
        try:
            data = data_ch1+1j*data_ch2               
        except:
            if type(data_ch1) is str:
                raise TypeError('data_ch1 is a string:\n',data_ch1)
            elif type(data_ch2) is str:
                raise TypeError('data_ch2 is a string:\n',data_ch2)
            else:
                raise CustomError("I don't know what's going on")
        if verbose: print 'checkpoint 2'
        return nddata(data,[-1],['t'],axis_coords = [self.time_axis])

    def retrieve_fast(self,verbose = False):
        #{{{ retrieve information about the points
        if verbose: 
            print 'Start time: \n'
            print self.time_start
        self.g.write('PAUS') #Stop data acquisition
        number_of_points = int(self.g.respond('SPTS?')) #Query Number of Points
        if verbose: print('Number of points: ' + str(number_of_points))
        frequency = int(round(2**(float(self.g.respond('SRAT?'))-4))) #Query frequency and convert to Hz
        if verbose: print('Frequency:'+str(frequency))
        #}}}
        #{{{ retrieve points for channels 1 and 2
        data_ch1 = self.g.respond('TRCL?' + '1,0,' + str(number_of_points),length = number_of_points*4,until = 'eoi',verbose = verbose) #Return Ch. 1
        data_ch2 = self.g.respond('TRCL?' + '2,0,' + str(number_of_points),length = number_of_points*4,until = 'eoi',verbose = verbose) #Return Ch. 2
        if verbose: print 'checkpoint 1'
        if verbose: print(len(data_ch1))
        #data_ch1 = data_ch1[::-1]
        #data_ch2 = data_ch2[::-1]
        #}}}
        #{{{ construct the structured data arrays and the time axis
        self.time_axis = r_[0:double(number_of_points)]/frequency+self.time_start
        if verbose: print 'Raw string:%s' %data_ch1
        #data_ch1 = fromstring(data_ch1,dtype = [('mantissa','<i2'),('exp','<u1'),('zero','<u1')])
        #data_ch2 = fromstring(data_ch2,dtype = [('mantissa','<i2'),('exp','<u1'),('zero','<u1')])
        data_ch1 = fromstring(data_ch1,dtype = [('mantissa','<i2'),('exp','<i2')])
        data_ch2 = fromstring(data_ch2,dtype = [('mantissa','<i2'),('exp','<i2')])
        #}}}
        #{{{ convert the data arrays
        if verbose: print 'Formatted string:';print data_ch1 
        data_ch1 = (array(data_ch1['mantissa'],dtype = 'float'))*(2**(array(data_ch1['exp'],dtype = 'float')-124.))
        if verbose: print 'Converted to Values:';print data_ch1
        data_ch2 = (array(data_ch2['mantissa'],dtype = 'float'))*(2**(array(data_ch2['exp'],dtype = 'float')-124.))
        #}}}
        #{{{ return the data as a complex nddata object
        try:
            data = data_ch1+1j*data_ch2               
        except:
            if type(data_ch1) is str:
                raise TypeError('data_ch1 is a string:\n',data_ch1)
            elif type(data_ch2) is str:
                raise TypeError('data_ch2 is a string:\n',data_ch2)
            else:
                raise CustomError("I don't know what's going on")
        if verbose: print 'checkpoint 2'
        return nddata(data,[-1],['t'],axis_coords = [self.time_axis])
        #}}}
    def start_storage(self):
        #self.front_display(1,'X') #These two lines were commented out for speed
        #self.front_display(2,'Y') #imporoved speed of function from 2.10599 sec to 0.002000 sec!
        self.g.write('REST')
        #self.time_start = time.time()
        self.g.write('STRT')
        self.time_start = time.time()
        return self.time_start
    def reset_buffer(self):
        self.g.write('REST')
#}}}
#{{{ ESP380 field controller
class field_controller ():
    def __init__(self,gpibtalk=4,gpiblisten=36,sock = 1234,ip = '127.0.0.1'):
        # talk means device talks, computer listens; listen means device listens, computer talks
        self.g = gpib(ip = ip,sock = sock)
        self.gpiblisten = gpiblisten
        self.gpibtalk = gpibtalk
        self.g.setaddr(self.gpiblisten)
        try:
            self.g.write('FC',addr = self.gpiblisten)
            stringfromfc = self.g.oldread(addr = self.gpibtalk)
            field = stringfromfc[22:-2]
        except:
            del self.g
            raise ValueError('Field controller is not responding!')
        if field == '':# try again
            self.g.write('FC',addr = self.gpiblisten)
            stringfromfc = self.g.oldread(addr = self.gpibtalk)
            field = stringfromfc[22:-2]
        if field == '':
            del self.g
            raise ValueError('Field controller is responding with the string "%s"!'%stringfromfc)
        else:
            print 'Connected to field controller. Current field:',field[3:],'G\n'
    def close(self):
        raise ValueError('I changed this --> you no longer call instance.close() and instead, you say del instance')
    def __del__(self):
        r'Close connection to field controller'
        del self.g
    def read_led(self):
        r'Read LEDs of field controller front display'
        self.g.write('LE',self.gpiblisten)
        return self.g.oldread(self.gpibtalk)[22:-2]
    def set_field(self,field):
        r'Set center field in Gauss'
        self.g.write('CF'+str(field),self.gpiblisten)
    def read_field(self):
        r'Read current field in Gauss'
        self.g.write('FC',self.gpiblisten)
        return self.g.oldread(self.gpibtalk)[22:-2]
    def read_sweep_address(self):
        r'Read current sweep address'
        self.g.write('SA',self.gpiblisten)
        return self.g.oldread(self.gpibtalk)
    def read_status_led(self):
        r'Read status LEDs'
        self.g.write('LE data',self.gpiblisten)
        return self.g.oldread(self.gpibtalk)
    def set_mode(self,mode):
        r'Set mode - mode must be integer 0, 2 or 3'
        self.g.write('MOD %s' %str(mode))
        return self.g.oldread(self.gpibtalk)
    def set_time(self,sweeptime):
        r'Set sweep time in seconds'
        self.g.write('TM'+str(sweeptime),self.gpiblisten)
    def set_width(self,sweepwidth):
        r'Set sweep width in Gauss'
        self.g.write('SW'+str(sweepwidth),self.gpiblisten)
    def sweep_up(self):
        r'Sweep upwards'
        self.g.write('SU',self.gpiblisten)
    def sweep_down(self):
        r'Sweep downwards'
        self.g.write('SD',self.gpiblisten)
    def stop(self):
        r'Stop sweep'
        self.g.write('ST',self.gpiblisten)
#}}}

#{{{ sampling scope
class sampling_scope ():
    def __del__(self):
       del self.g
    def __init__(self,gpibaddress=1,sock = 1234,ip = '127.0.0.1'):
        self.g = gpib(ip = ip)
        self.gpibaddress = gpibaddress
        self.g.setaddr(self.gpibaddress)
        print self.id()
        return
    def capture(self,
            do_plot = True,
            input_center_frequency = 9.6e9,
            bandpass2 = 1e9,
            bandpass = 5e9,
            detbyphase = False, #determine frequency by summing up the phases, rather than picking the peak of the FT
            channel = 2,
            navg = 64
            ):
        readwaveform_kwargs = {}
        if navg is not None:
            readwaveform_kwargs.update({'navg':navg})
        self.g.write('NAVG %d'%(navg/3))
        print "How many averages (/3) ?:",self.g.respond('NAVG?')
        #{{{ loop to acquire set number of averages
        self.g.write('CONDacq TYPE:AVG')
        n = 100000
        while n>0:
            n = int(self.g.respond('CONDacq? REMAINING').split(':')[1])
            #print "I have %d averages left"%n
        #}}}
        #{{{ now, literally just read the waveforms
        plain_yig_input = self.read_waveform(channel=1) # grab plain yig
        signaldata_input = self.read_waveform(channel=channel) # grab waveform
        #}}}
        plain_yig = plain_yig_input.copy()
        data = signaldata_input.copy()
        fl = figlist() # new figure list
        ### data = self.yig_as(signaldata, plain_yig) # create analytic signal ### OLD METHOD ###
        #{{{ apply lowpass and convert to analytic signal in one shot
        plain_yig.ft('t', shift=True)
        data.ft('t', shift=True)
        #{{{ make a new analytic signal array that I use for frequency determination
        #{{{ show what the FT looks like
        if do_plot:
            fl.next('yigft')
            plot(abs(plain_yig),label='plain yig')
            plot(abs(data),label='iq output')
            xlabel('f / GHz')
            legend()
        #}}}
        #plain_yig['t',lambda x: x<0] = 0.
        center_frequency = abs(plain_yig['t',lambda x: x>0]).argmax('t').data
        #print('the peak of the YIG FT is at %0.2f GHz'%(center_frequency/1e9))
        #}}}
        if center_frequency < input_center_frequency - bandpass:
            print "WARNING (capture): You tell me the frequency is %0.2f GHz, but I have determined it to be %0.2f GHz, so I'm assuming that your data is aliased, so I'm pulling the negative peak rather than the positive peak for the analytic signal"%(input_center_frequency/1e9,center_frequency/1e9)
            center_frequency *= -1
        plain_yig['t',lambda x:logical_or(x>bandpass+center_frequency,x<center_frequency-bandpass)] = 0.
        data['t',lambda x:logical_or(x>bandpass+center_frequency,x<center_frequency-bandpass)] = 0.
        plain_yig.ift('t', shift=True)
        data.ift('t', shift=True)
        #}}}
        #{{{ convert phase of test analytic signal to difference from reference
        data.data = data.data.copy() / plain_yig.data * abs(plain_yig.data) # subtract the phase of the reference waveform from the test waveform
        data.ft('t', shift = True)
        data['t',lambda x: abs(x)>bandpass2] = 0.
        data.ift('t', shift = True)
        #}}}
        if do_plot:
            data_plot = data.copy()
            fl.next('as') # start on a plot called ``analytic signal''
            #x = data_plot.getaxis('t')
            #x /= 1e-9
            #data_plot.rename('t','t / ns')
            # plot(data_plot,alpha = 0.1,label = 'original')
            plot(abs(data_plot),linewidth = 2, label = 'amplitude')
            ax1 = gca() # amplitude axis
            # ax1.set_ylim([-0.015,0.015])
            twinx() # twin y axes
            ax2 = gca() # phase axis
            angles = data_plot.copy() # generate phase data
            angles.data = angle(angles.data)*180./pi
            plot(angles,'r.',markersize = 2,ax = ax2, label = 'phase')
            ax1.set_ylim(0)
            ax2.set_ylim([-180,180])
            autolegend(ax = ax1,ax2 = ax2)
        return data
    def yig_as(self, input_data, input_plainyig,
            maxdispfreq = 20e9,
            lowpass = 3e9,
            centerfreq = None,
            filterwidth = 3.0e9,
            threshold = 0.05,
            verbose = False,
            verbose2 = True,
            center_time = 1.8e-9): # return mixed down waveform with absolute phase information
        #{{{ analyze plain yig
        data = input_plainyig.copy()
        data.ft('t',shift = True)
        x = data.getaxis('t') # to reference x axis
        data.data[x<0.0] = 0 # make analytic signal by dropping negative frequencies
        normdata = abs(data)
        normdata.data[abs(x)<lowpass] = 0 # apply lowpass, default freq: 3e9
        normdata.data[normdata.data<threshold*max(normdata.data)] = 0 # cut off everything that is less than x% of maximum, default threshold: 5%
        #{{{ normalize abs of data for frequency determination
        avgdata = normdata.copy().mean('t').data # average magnitude of the ft
        normdata /= avgdata # normalize to 1
        #}}}
        # determine mean frequency, then ift and mix down
        meanfreq = mean(normdata.data*x)
        if verbose2:
            print ' '
            print 'meanfreq = %0.4f GHz'%(meanfreq/1e9)
        data.ift('t',shift = True) # inverse ft
        x = data.getaxis('t')
        data.data *= exp(-1j*2*pi*meanfreq*x) # mix down
        # find angles
        angles = data.copy()
        angles.data /= abs(angles.data) # normalize
        #{{{ avoid rollover
        for i in arange(len(angles.data)-1):
            angles.data[i] = angles.data[i+1] / angles.data[i] # get e^i(phi_2-phi_1)
        angles.data = angle(angles.data) # generate phase difference data
        if verbose: print angles.data[1], 'new 1'
        angles.data = cumsum(angles.data) # cumulative sum to get phase data back
        yangles = angles.data[100:-100] # y-axis, drop first and last 100 items for precision
        xangles = angles.getaxis('t')[100:-100] # same for x-axis
        #}}}
        (corrfreq,b1) = polyfit(xangles,yangles,1) # linear regression
        corrfreq /= 2*pi
        if verbose2:
            print 'correction frequency = %0.4f MHz'%(corrfreq/1e6) # print slope = correction frequency
        meanfreq += corrfreq
        if verbose2:
            print 'new meanfreq = %0.4f GHz'%(meanfreq/1e9)
        data.data *= exp(-1j*2*pi*corrfreq*x)
        angles.data = angle(data.data)*180./pi
        phasecorr = mean(angles.data)/180.*pi # determine phase correction
        #}}} analyze plain yig
        #{{{ analytic signal
        data = input_data.copy()
        data.ft('t',shift = True)
        x = data.getaxis('t') # to reference x axis
        data.data[x<0.0] = 0 # make analytic signal by dropping negative frequencies
        data.data[abs(x-meanfreq)>filterwidth] = 0 # filter
        data.ift('t',shift = True) # inverse ft
        x = data.getaxis('t') # to reference x axis
        data.data *= exp(-1j*(2*pi*meanfreq*x + phasecorr)) # generate complex signal with right frequency and phase correction
        return data
    def clear(self):
        string = self.g.write('CLEar ALLTrace')
    def command(self,cmd):
        if cmd[-1] == '?':
            reply = self.g.respond(cmd)
            return reply
        else:
            self.g.write(cmd)
    def id(self):
        try:
            idstring = self.g.respond('ID?')
        except:
            raise ValueError('The sampling scope is not responding to the ID command --> maybe it\'s not hooked up, or you\'re talking to the wrong GPIB address\n\n OR you LIKELY need to set the terminator to EOI/LF --> really, I should figure out how to program this')
        return 'Connected to '+idstring
    def init_active_cancel(self):
        # initializes sampling scope with the settings we need for active cancellation
        self.timebase(20e-9)
        self.resolution(5120)
        self.position(516.9e-9)
        #self.g.write('REM TRA1')
        #self.g.write('REM TRA2')
        self.g.write('TRA1 DES:"AVG(M1)"')
        self.g.write('NAVG 256')
        self.g.write('TRA2 DES:"AVG(M2)"')
        self.g.write('NAVG 256')
    def position(self, pos):
        # defines position after trigger pulse in seconds
        self.g.write('MAINPos %0.6e'%pos)
    def read_waveform(self, channel = 1, verbose = False):
        # reads out waveform of specified channel and returns nddata
        junk = self.g.respond('ID?') # to clear any junk
        self.g.write('OUTPut TRAce' + str(channel))
        mystring = self.g.respond('WAVFRM?')
        if verbose: print 'last element of string is',repr(mystring[-1])
        if verbose: print 'string is',mystring
        trace_re = re.compile(r'CURVE CRVID:TRACE([0-9]+),(.+)')
        header_re = re.compile(r'WFMPRE WFID:TRACE([0-9]+),(.+)')
        found = True
        for thisstring in mystring.split(';'):
            m = header_re.match(thisstring)
            if m:
                index = int(m.groups()[0])
                if verbose: print 'index of header',index
                vars = map(lambda x: tuple(x.split(':',1)), m.groups()[1].split(','))
                vars = dict(vars)
            m = trace_re.match(thisstring)
            if m:
                index = int(m.groups()[0])
                if verbose: print 'index of waveform',index
                data = array(map(double,m.groups()[1].split(',')))
        returndata = nddata(data*double(vars['YMULT']),[-1],['t'],axis_coords = [double(vars['XINCR'])*r_[0:size(data)]])
        returndata.set_units('t','s')
        return returndata
    def resolution(self, datapoints):
        # sets number of datapoints (512 to 5120)
        if type(datapoints) is str:
            if datapoints == 'max':
                self.g.write('TBMain Length:5120')
            elif datapoints == 'min':
                self.g.write('TBMain Length:512')
            else:
                raise ValueError('Only accepted strings are max and min.')
        else:
            self.g.write('TBMain Length:%i'%datapoints)
    def timebase(self, resolution):
        # sets timebase in s/div
        self.g.write('TBMain Time:%0.1e'%resolution)
>>>>>>> public
#}}}
