class nddata_hdf5(nddata):
    def __repr__(self):
        if hasattr(self, "_node_children"):
            return repr(self.datanode)
        else:
            return nddata.__repr__(self)
        atexit.register(self._cleanup)

    def _cleanup(self):
        if hasattr(self, "_node_children"):
            self.h5file.close()
            del self.h5file
            del self.datanode
        return

    def __init__(self, pathstring, directory="."):
        self.pathstring = pathstring
        # try:
        self.h5file, self.datanode = h5nodebypath(
            pathstring, check_only=True, directory=directory
        )
        logger.debug("about to call _init_datanode")
        self._init_datanode(self.datanode)
        atexit.register(self._cleanup)

    def _init_datanode(self, datanode, **kwargs):
        datadict = h5loaddict(datanode)
        # {{{ load the data, and pop it from datadict
        try:
            datarecordarray = datadict["data"][
                "data"
            ]  # the table is called data, and the data of the table is called data
            mydata = datarecordarray["data"]
        except Exception:
            raise ValueError("I can't find the nddata.data")
        try:
            kwargs.update({"data_error": datarecordarray["error"]})
        except Exception:
            logger.debug(strm("No error found\n\n"))
        datadict.pop("data")
        # }}}
        # {{{ be sure to load the dimlabels
        mydimlabels = [j.decode("utf-8") for j in datadict["dimlabels"]]
        if len(mydimlabels) == 1:
            if len(mydimlabels[0]) == 1:
                mydimlabels = list(
                    [mydimlabels[0][0]]
                )  # for some reason, think I need to do this for length 1
        # }}}
        # {{{ load the axes and pop them from datadict
        datadict.pop("dimlabels")
        if "axes" in list(datadict.keys()):
            myaxiscoords = [None] * len(mydimlabels)
            myaxiscoordserror = [None] * len(mydimlabels)
            myaxis_units = [None] * len(mydimlabels)
            logger.debug(
                strm(
                    "about to read out the various axes:",
                    list(datadict["axes"].keys()),
                )
            )
            for axisname in list(datadict["axes"].keys()):
                try:
                    axisnumber = mydimlabels.index(axisname)
                except AttributeError:
                    raise AttributeError(
                        strm(
                            "mydimlabels is not in the right format!\nit looks like this:\n",
                            mydimlabels,
                            type(mydimlabels),
                        )
                    )
                except ValueError:
                    raise ValueError(
                        strm(
                            "mydimlabels is not in the right format!\nit looks like this:\n",
                            mydimlabels,
                            type(mydimlabels),
                        )
                    )
                recordarrayofaxis = datadict["axes"][axisname]["data"]
                if "axis_coords_units" in datadict["axes"][axisname].keys():
                    myaxis_units[axisnumber] = datadict["axes"][axisname][
                        "axis_coords_units"
                    ]
                else:
                    if ("Scans" in axisname) or ("ph" in axisname):
                        pass
                    else:
                        print(
                            "You didn't set units for %s before saving the data!!!"
                            % axisname
                        )
                myaxiscoords[axisnumber] = recordarrayofaxis["data"]
                if "error" in recordarrayofaxis.dtype.names:
                    myaxiscoordserror[axisnumber] = recordarrayofaxis["error"]
                datadict["axes"][axisname].pop("data")
                for k in list(datadict["axes"][axisname].keys()):
                    logger.debug(
                        strm(
                            "Warning, attribute",
                            k,
                            "of axis table",
                            axisname,
                            "remains, but the code to load this is not yet supported",
                        )
                    )
                datadict["axes"].pop(axisname)
            kwargs.update({"axis_coords": myaxiscoords})
            kwargs.update({"axis_coords_units": myaxis_units})
            kwargs.update({"axis_coords_error": myaxiscoordserror})
        elif len(mydimlabels) > 1:
            raise ValueError(
                "The current version uses the axis labels to"
                "figure out the shape of the data\nBecause you stored"
                "unlabeled data, I can't figure out the shape of the"
                "data!!"
            )
            # the reshaping this refers to is done below
        # }}}
        logger.debug(
            strm(
                "about to initialize data with shape",
                mydata.shape,
                "labels",
                mydimlabels,
                "and kwargs",
                kwargs,
            )
        )
        nddata.__init__(self, mydata, mydata.shape, mydimlabels, **kwargs)
        # {{{ reshape multidimensional data to match the axes
        if len(mydimlabels) > 1:
            det_shape = []
            for thisdimlabel in mydimlabels:
                try:
                    temp = self.getaxis(thisdimlabel)
                except Exception:
                    temp = -1  # no axis is given
                if isinstance(temp, np.ndarray):
                    temp = len(temp)
                det_shape.append(temp)
            try:
                self.data = self.data.reshape(det_shape)
            except Exception:
                raise RuntimeError(
                    strm(
                        "The data is of shape",
                        self.data.shape,
                        "and I try to reshape it into",
                        det_shape,
                        "corresponding to the dimensions",
                        mydimlabels,
                        "--> this fails!",
                    )
                )
        # }}}
        for remainingattribute in list(datadict.keys()):
            self.__setattr__(remainingattribute, datadict[remainingattribute])
        self.h5file.close()
        del self.h5file
        del self.datanode
        return


