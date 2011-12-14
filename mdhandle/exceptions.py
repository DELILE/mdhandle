# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Contains custom exceptions for mdhandle module.

.. note::
   v0.9 not implemented.

"""

#TODO: Integration with logging functionality.
#TODO: Replace existing general Exception and UserWarning calls with exceptions
#      or logging.



#------------------------------------------------------------------------------#
#    Exceptions
#------------------------------------------------------------------------------#

class Error(Exception):
    """
    General mdhandle exception.
    """
    pass


class InputError(Error):
    """
    Raised when input parameters are wrong format, size, type
    """
    def __init__(self, *args):
        self.args = args
