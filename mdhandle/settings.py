# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Contains basic settings (system, LAMMPS format, file format,etc.) for
:mod:`mdhandle` package.

"""
# TODO: Default file name format?
    # see utilities.time_from_filename(...) for info
    # REGEX_time = regex definition or function to define

import sys
import numpy as np

# ----------------------------------------------------------------------------
#       SITE SETTINGS
# ----------------------------------------------------------------------------

# Site-wide settings can be modified here, or at runtime.

# specifying program to use for .gz files
GUNZIP = 'gunzip'

# Location for temporary files.
TMP_LOCATION = '/tmp'

# Structure of LAMMPS ASCII dump header.  One line per list entry.
LAMMPSHEADER_LEGACY = ['ITEM: TIMESTEP', 'INT_TIME', 'ITEM: NUMBER OF ATOMS',
                  'INT_NUM_ATOMS', 'ITEM: BOX BOUNDS', 'FLOAT_XLO FLOAT_XHI',
                  'FLOAT_YLO FLOAT_YHI', 'FLOAT_ZLO FLOAT_ZHI', 'ITEM: ATOMS']

LAMMPSHEADER = ['ITEM: TIMESTEP', 'INT_TIME', 'ITEM: NUMBER OF ATOMS',
                    'INT_NUM_ATOMS', 'ITEM: BOX BOUNDS', 'FLOAT_XLO FLOAT_XHI',
                    'FLOAT_YLO FLOAT_YHI', 'FLOAT_ZLO FLOAT_ZHI',
                    'ITEM: ATOMS <DUMP-OUTPUT-COLUMNS>']

# Defining the options available for atomwise or gridded data.
PROP_SIZE = {'scalars': 1, 'vectors': 3, 'tensors': 9, 'symm_tensors': 6}
DTYPES = ['int', 'float']
LOCATIONS = ['Node', 'Cell']

# Required columns within any column mapping (mdhandle.readers.column_mapping)
REQUIRED_DATA = np.array(['id', 'x', 'y', 'z'])

# Name of default atomwise dataset within a snapshot
DEFAULT_DATASET = 'rawSimResults'

#=============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
