# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Provides utility functions for manipulating tensor collections. In this
context, tensor collections are (Nx6) or (Nx9) :class:`numpy.arrays` where each row is a tensor associated with an individual atom.

The order of components in a 6-component symmmetric tensors is:
    ``00, 01, 02, 11, 12, 22``

    ``xx, xy, xz, yy, yz, zz``

The order of components in a 9-component general tensor is:
    ``00, 01, 02, 10, 11, 12, 20 21, 22``

    ``xx, xy, xz, yx, yy, yz, zx, zy, zz``

"""

# TODO: Vectors along surfaces
# TODO: Principle directions, eigenvalues (invariants)
# TODO: Consider refactoring so that each atomwise item a 2D matrix object
#       rather than a row vector.

#------------------------------------------------------------------------------

import sys

import numpy as np

from mdhandle.logger import Logger

#------------------------------------------------------------------------------

logger = Logger()

#------------------------------------------------------------------------------


def grow(cont_tens):
    """
    Returns a tensor collection grown from input list/tuple, or comma separated
    list of existing vector collections.

    Parameters
    -----------
    cont_tens : iterable of :class:`numpy.ndarray`, (Nx6 or Nx9)
        List or tuple of tensor collections to be stacked.

    Returns
    --------
    grown : :class:`numpy.ndarray`, ((N1 + N2 + ...)x6) or ((N1 + N2 + ...)x9)
        Returns new tensor collection.

    See Also
    ---------
    :meth:`mdhandle.vectors.grow` for related example.

    """
    for i in cont_tens:
        assert isinstance(i, np.ndarray)
        if i.ndim == 1:
            assert (i.shape == (6,) or i.shape == (9,))
        elif i.ndim == 2:
            assert (i.shape[1] == 6 or i.shape[1] == 9)
    grown = np.vstack(tuple(cont_tens))
    return grown


def build_from_scalars(cont_tens):
    """
    Returns a tensor collection ((Nx6) or (Nx9) :class:`numpy.ndarray`) 
    from its scalar components ((Nx1) :class:`numpy.ndarray`).

    Parameters
    ----------
    cont_vecs : iterable of :class:`numpy.ndarray`, (Nx1)
        Collection of constituent scalars: ``(t0, t1, t2, t3, t4, t5, ...)``.
        Each scalar contains atomwise scalar.

    Returns
    --------
    new_tens : :class:`numpy.ndarray`, (Nx6) or (Nx9)
        Returns tensor collection in the form
        ``numpy.array([t0, t1, t2, t3, t4, t5, ... ])``.

    See Also
    ---------
    :meth:`mdhandle.vectors.build_from_scalars` for related example.

    """
    assert (len(cont_tens) == 6 or len(cont_tens) == 9)

    shape_test = cont_tens[0].shape
    for t_item in cont_tens:
        assert isinstance(t_item, np.ndarray)
        # ASSERT: All scalar components are the same length
        assert t_item.shape == shape_test

    new_tens = np.column_stack(cont_tens)
    return new_tens


def is_symmetric(tens):
    """
    Tests for the symmetry of the supplied tensor.

    Parameters
    ----------
    tens : ``numpy.ndarray``   (Nx6 or Nx9)
        Tensor collection.

    Returns
    -------
    return_value : boolean
        If ``True``, ``tens`` is symmetric.

    """
    return_value = False
    
    if tens.ndim == 1:
        if len(tens) == 6:
            return_value =  True
        elif len(tens) == 9:
            if ((tens[1] == tens[3]) and (tens[2] == tens[6])
                    and (tens[5] == tens[7])):
                return_value = True
            else:
                return_value = False
        else:
            logger.warning('Tensor is improper shape: %s' % tens.shape)
            return_value = False

    elif tens.ndim == 2:
        if tens.shape[0,:] == 6:
            # 6-component tensor collections are symmetry by default.
            return_value = True
        elif tens.shape[0,:] == 9:
            if ( np.allclose(tens[:,1], tens[:,3])
                            and np.allclose(tens[:,2], tens[:,6])
                            and np.allclose(tens[:,5], tens[:,7]) ):
                logger.user_message('9-component tensor is symmetric.')
                return_value = True
            else:
                return_value = False
        else:
            logger.warning('Tensor is improper shape. Shape: %s' % tens.shape)
            return_value = False
    else:
        logger.user('Tensor is >2 dimensional array. Dims: %s' % tens.ndim)
        return_value = False
    
    return return_value


#=============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
