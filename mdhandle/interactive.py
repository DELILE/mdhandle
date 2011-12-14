# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Enables ipython interaction for from automated scripts and other routines.

Calling :func:`interactive.launch_ipshell` will launch iPython as a subshell.


**References**:

* http://ipython.org/ipython-doc/dev/interactive/reference.html#embedding-ipython

"""

#----------------------------------------------------------------------------

import sys

# Import the embed function
from IPython.frontend.terminal.embed import InteractiveShellEmbed

# -----------------------------------------------------------------------------

def launch_ipython():
    """
    Launches iPython shell, or sub-shell if already in iPython.

    Wrapper around ``InteractiveShellEmbed`` function from
    ``IPython.frontend.terminal.embed``.

    """
    # Trying if already in iPython 
    try:
        get_ipython
        # TODO: returning None required to account for bug as in:
        #       https://github.com/ipython/ipython/issues/1045
        #       Cannot embed IPython in IPython.
        return None
    except NameError:
        banner = exit_msg = ''
    else:
        banner = '\n*** Entering Nested interpreter ***\n'
        exit_msg = '\n*** Back in main IPython ***\n'

    ipshell = InteractiveShellEmbed(banner1=banner, exit_msg=exit_msg)
    return ipshell()

# ============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
