# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Helper functions for working with atom and gridded coordinates as of [Nov 2012]

TODO: 
    - Add conservative or not flag to convert_to_grid_idx()
        + This would take care of +/- 1 to deal with end point off by one
          issues due to Numpy slicing.

"""
import numpy as np
import glob

from mdhandle import snap


def convert_to_grid_idx(input_value, grid_spacing, direction, zero_value=2.):
    """
    This function converts coordinates from atom coordinates
    to gridded coordinates.
    
    zero_value: Minimum value in the gridded data in the direction.
    
    """
    zero_value = float(zero_value)
    
    if direction == 'x':
        result = np.floor((float(input_value) - 6. - zero_value) / grid_spacing)

    if direction == 'y':
        result = np.floor( (float(input_value) - zero_value) / grid_spacing)

    if direction == 'z':
        result = np.floor((float(input_value) + 6. - zero_value) / grid_spacing)

    return int(result)
    
def convert_from_grid_idx(input_idx, grid_spacing, direction, zero_value=2.):
    """
    zero_value: Minimum value in the gridded data in the direction.
    
    """
    zero_value = float(zero_value)
    
    if direction == 'x':
        result = float(input_idx*grid_spacing + 6. + zero_value)
    if direction == 'y':
        result = float(input_idx*grid_spacing + zero_value)
    if direction == 'z':
        result = float(input_idx*grid_spacing - 6. + zero_value)
    return int(result) 
    

def grid_correction(input):
    return input+np.array([6., 0., -6.])

def flip_vectors(indep,dependent):
    """
    Flips vectors end to end for plotting purposes.
    
    Use: Flipping output from matplotlib contour line generation.
    
    """
    if np.any(np.diff(indep)[0] < 0):
        return indep[-1::-1], dependent[-1::-1]
    else:
        return indep, dependent
        
# -----------------------------------------------------------------------------

def timestep_num_from_frame_number(num):
    flist = glob.glob('*.h5')
    
    s = snap.Snap(flist[num])
    s.gather_data()
    return s.meta['time']
    
def timestep_from_nd_time(nd_time):
    flist = glob.glob('*.h5')
    
    print('Assuming constant timestep size')
    
    s = snap.Snap(flist[0])
    s.gather_data()
    return int(nd_time / s.meta['timestep'])    