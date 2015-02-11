# note that there were two different versions of this code -- the one a few commits before the comment is introduced is the "older" version that was used in CNSI
import os
#import termios
import socket
import time
from struct import *
from pylab import *
from socket import socket
import scipy.io as mio
DEFAULTCOM = 4

#{{{ the underlying gpib class, which has older instruments integrated (separate them later
class gpib():
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

#{{{ this section gives wrappers for specific instruments
class eip_powermeter ():
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
        self.g.write(self.gpibaddress,'DA')# turn on the display
        self.g.write(self.gpibaddress,'HP')# if we don't do this, the display freezes
        self.g.write(self.gpibaddress,'RP')# no longer output only mode
        self.g.write(self.gpibaddress,'FP')# turn off "fast mode"??
        self.g.write(self.gpibaddress,'R0')# switch back to high res
        self.g.close()
class hp8672a():
    def __init__(self,comport = DEFAULTCOM,gpibaddress = 3):
        self.g = gpib(comport - 1) # on port 4 with the current connector
        self.gpibaddress = gpibaddress
        self.stepsize = 0.5e6 # this is a lie, but it's used as a default by wobbandmin
    def set_frequency(self,frequency):
        self.g.write(self.gpibaddress,'P%08dZ0'%int(round(frequency*1e-3)))# just use the 10 GHz setting, and fill out all the other decimal places with zeros
        return
    def close(self):
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
        #self.g.write(self.gpibaddress,'DE')         # Display Enable
     
        
##        self.g.write(self.gpibaddress,'HP')# if we don't do this, the display freezes
##        self.g.write(self.gpibaddress,'RP')# no longer output only mode
##        self.g.write(self.gpibaddress,'FP')# turn off "fast mode"??
##        self.g.write(self.gpibaddress,'R0')# switch back to high res
        self.g.close()

#}}}
