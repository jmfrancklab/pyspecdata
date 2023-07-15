import numpy as np

class CustomError(Exception):
    def __init__(self, *value, **kwargs):
        raise NotImplementedError("You should get rid of CustomError and use explain_error instead")
        return
def emptytest(x): # test is it is one of various forms of empty
    if x is None: return True
    if np.size(x) == 0: return True
    if np.size(x) == 1 and x == None: return True
    return False
def explain_error(e):
    '''Allows you to wrap existing errors with more explanation

    For example:

    >    except BaseException  as e:
    >        raise IndexError("I can't find the node "+pathstring+explain_error(e))
    >                + '\n'.join(['>\t'+j for j in str(e).split('\n')]))# this indents
    '''
    exc_type,exc_obj,exc_tb = exc_info()
    #code_loc = strm(os.path.relpath(exc_tb.tb_frame.f_code.co_filename,os.getcwd()), 'line', exc_tb.tb_lineno)
    code_loc = strm(exc_tb.tb_frame.f_code.co_filename, 'line', exc_tb.tb_lineno)
    return ('\n> Original error (%s -- %s):\n'%(exc_type,code_loc)
            + '\n'.join(['>\t'+j
                for j in str(e).split('\n')]))# this indents
def maprep(*mylist):
    mylist = list(mylist)
    for j in range(0,len(mylist)):
        if not isinstance(mylist[j], str):
            mylist[j] = mylist[j].__repr__()
    return ' '.join(mylist)
