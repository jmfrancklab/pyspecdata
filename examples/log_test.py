"""
Logging Functionality
=====================

a demonstration of the logging functionality -- output is logged to
`~/pyspecdata.log` (or the same name with a number included if multiple scripts
are running at the same time"""

from pyspecdata import *
mylogger = init_logging("info")
mylogger.info("Something, something, something, dark side...")

