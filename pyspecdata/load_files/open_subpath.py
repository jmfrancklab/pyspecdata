from ..core import *
from ..datadir import dirformat
import os.path
from zipfile import ZipFile
import codecs

utf8reader = codecs.getreader('utf-8')

def open_subpath(file_reference,*subpath,**kwargs):
    """
    Parameters
    ----------
    file_reference: str or tuple
        If a string, then it's the name of a directory.
        If it's a tuple, then, it has three elements: the ZipFile object, the
        filename of the zip file (for reference), and the name of the file we're interested in
        within the zip file.
    test_only: bool
        just test if the path exists
    """
    logger.debug(strm("trying to open subpath: ",*subpath,kwargs))
    mode,test_only = process_kwargs([('mode','r'),
        ('test_only',False),
        ],kwargs)
    if isinstance(file_reference,str):
        logger.debug("this appears to be a standard (non-zip) file")
        if test_only:
            full_path = os.path.join(file_reference, *subpath)
            if os.path.exists(full_path):
                return True
            else:
                return False
        else:
            fp = open(os.path.join(file_reference,*subpath),mode)
    else:
        if isinstance(file_reference, tuple):
            logger.debug("this appears to be a zip file")
            if len(file_reference) == 3 and isinstance(file_reference[0], ZipFile):
                zf = file_reference[0]
                zip_basename = file_reference[1]
                name_inside_zip = file_reference[2]
                subfile = '/'.join((name_inside_zip,)+subpath)
                if test_only:
                    if subfile in zf.namelist():
                        return True
                    else:
                        return False
                if subfile in zf.namelist():
                    if 'b' in mode:
                        return zf.open(subfile)
                    else:
                        return utf8reader(zf.open(subfile))
                else:
                    raise ValueError(subfile+" not found in zip file")
            else:
                raise ValueError("open_subpath doesn't understand the format of the tuple passe to file_reference")
        else:
            raise ValueError("open_subpath doesn't understand the type of the file_reference")
    return fp
