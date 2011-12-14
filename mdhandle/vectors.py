# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Provide utility functions for manipulating vector collections.
In this context, vector collections are Nx3 :class:`numpy.ndarray` where N is the total number of atoms.  

Each row in the vector collection is a vector associated with the individual atom.

"""

import sys

import numpy as np

from mdhandle.logger import Logger

#------------------------------------------------------------------------------

logger = Logger()

#------------------------------------------------------------------------------


def grow(cont_vecs):
    """
    Returns a vector collection grown from input list/tuple, or comma separated
    list of existing vector collections.

    Parameters
    -----------
    cont_vecs : iterable of :class:`numpy.ndarray`, (Nx3)
        List or tuple of vector collections (:class:`numpy.ndarray`, (Nx3)) to 
        be stacked.

    Returns
    --------
    grown : :class:`numpy.ndarray` ((N1 + N2 + ...)x3 )
        Returns new vector collection (numpy.ndarray).

    Examples
    ---------
    >>> a = np.array([[   2.25313,   22.4553 ,   11.1074 ],
                      [ 149.999  ,    4.65435,   21.8793 ]])
    >>> b = np.array([ 153.184  ,    7.48026,  143.027  ])

    >>> grow((a,b))
    array([[   2.25313,   22.4553 ,   11.1074 ],
           [ 149.999  ,    4.65435,   21.8793 ],
           [ 153.184  ,    7.48026,  143.027  ]])

    """
    for i in cont_vecs:
        assert isinstance(i, np.ndarray)
        if i.ndim == 1:
            assert i.shape == (3,)
        elif i.ndim == 2:
            assert i.shape[1] == 3
    grown = np.vstack(tuple(cont_vecs))
    return grown


def build_from_scalars(cont_vecs):
    """
    Returns a vector collection (:class:`numpy.ndarray`, (Nx3)) from its scalar
    components (:class:`numpy.ndarray`, (Nx1)).

    Parameters
    ----------
    cont_vecs : iterable of :class:`numpy.ndarray`, (Nx1)
        Collection of constituent scalars: ``(v0, v1, v2)``.
        Each scalar contains atomwise scalar value.

    Returns
    --------
    vec_coll : :class:`numpy.ndarray`, (Nx3)
        Returns vector collection in the form
        ``numpy.array([v_x, v_y, v_z])``.

    Examples
    ---------
    >>> v0 = np.array([   2.25313,  149.999  ,   34.7775 ,  153.184  ])
    >>> v1 = np.array([ 22.4553 ,   4.65435,  59.7746 ,   7.48026])
    >>> v2 = np.array([  11.1074,   21.8793,   82.1284,  143.027 ])
    >>> build_from_scalars((v0, v1, v2))
    array([[   2.25313,   22.4553 ,   11.1074 ],
           [ 149.999  ,    4.65435,   21.8793 ],
           [  34.7775 ,   59.7746 ,   82.1284 ],
           [  153.184 ,   7.48026  ,  143.027 ]])

    """

    # Will raise exception if cannot be unpacked
    (v0, v1, v2) = cont_vecs

    assert isinstance(v0, np.ndarray)
    assert isinstance(v1, np.ndarray)
    assert isinstance(v2, np.ndarray)

    # ASSERT: All scalar components are the same length
    assert v0.shape == v1.shape
    assert v0.shape == v2.shape
    assert v1.shape == v2.shape

    vec_coll = np.column_stack(cont_vecs)
    return vec_coll


def magnitude(vec):
    """
    Returns magnitude of each row within a vector collection.

    Parameters
    ----------
    vec : :class:`numpy.ndarray`, (Nx3)
        Stack of vectors, a vector collection.

    Returns
    -------
    mags : :class:`numpy.ndarray`, (Nx1)
        Magnitude of each row vector within the vector collection.

    Examples
    --------
    >>> vec = np.array([[ 1 , 2 , 4],
                        [ 2 , 3 , 6],
                      ])
    >>> magnitude(vec)
    array([4.5825, 7. ])

    """
    if vec.ndim == 1:
        mags = np.linalg.norm(vec)
    elif vec.ndim == 2:
        mags = np.sqrt(sq_magnitude(vec))
    return mags


def sq_magnitude(vec):
    """
    Returns the square magnitude of each row within a vector collection.

    Parameters
    ----------
    vec : :class:`numpy.ndarray`, (Nx3)
        Stack of vectors, a vector collection.

    Returns
    --------
    sq_mag : :class:`numpy.ndarray` (Nx1)
       Square magnitude of each row vector within the vector collection.

    Examples
    --------
    >>> vec = np.array([[ 1 , 2 , 4],
                        [ 2 , 3 , 6],
                      ])
    >>> sq_magnitude(vec)
    array([ 21.,  49.])

    """
    if vec.ndim == 1:
        sq_mag = np.power(np.linalg.norm(vec), 2)
    elif vec.ndim == 2:
        sq_mag = (vec*vec).sum(axis=1)
    return sq_mag


def unit_vec(vec):
    """
    Returns collection of unit vectors derived from input vector collection.

    Parameters
    -----------
    vec : :class:`numpy.ndarray` (Nx3)
        Stack of vectors, a vector collection.

    Returns
    --------
    unit_v : :class:`numpy.ndarray`, (Nx3)
        Collection of unit vectors.

    Examples
    --------
    >>> vec = array([[ 1 , 2 , 4],
                     [ 2 , 3 , 6],
                   ])
    >>> unit_vec(vec)
    array([[ 0.218..., 0.436..., 0.872...],
           [ 0.286..., 0.428..., 0.857...],
         ])

    """
    mag = magnitude(vec)

    if vec.ndim == 1:
        unit_v = vec / mag
    elif vec.ndim == 2:
        unit_v = vec / build_from_scalars((mag, mag, mag))
    return unit_v


def vector_type_conv(vec, new_dtype):
    """
    Converts vector collection into required data format.
    Currently converts only to: ``float``, ``int``

    Parameters
    ----------
    vec : :class:`numpy.ndarray`, (Nx3)
        Vector collection.
    new_dtype : string, {'int' | 'float' }
        Desired type of data.

    Returns
    --------
    vec : :class:`numpy.ndarray` (Nx3)
        ``vec`` with modified dtype.

    """
    new_dtype = new_dtype.lower()
    # TODO: Move list of acceptable dtypes to settings.  System dependent.
    if new_dtype not in ['int',   'float',
                         'int64', 'float64',
                         'int32', 'float32']:
        raise UserWarning("Desired data type not supported.  \
                                                    Must be 'int', 'float'.")

    if vec.dtype != np.dtype(new_dtype):
        vec = vec.astype(new_dtype)
    # If dtype is already ok - do nothing and return original vec
    else:
        pass

    return vec

# =============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
