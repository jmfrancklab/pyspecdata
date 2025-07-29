#ifndef PyDataType_SET_ELSIZE
#  define PyDataType_SET_ELSIZE(descr, size) \
          ((descr)->elsize = (npy_intp)(size))
#endif
