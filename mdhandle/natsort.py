# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Provides functions for natural string sorting (ie. strings with letters
and numbers).

**Credits**:

* Seo Sanghyeon
* Connelly Barnes

"""

# TODO: Manipulations could be done more smoothly with list list comprehensions
# in py >=2.5.

import re
import copy

def try_int(s):
    """
    Internal function to :mod:`natsort`.  Called by :func:`natsort_key`.

    Parameters
    -----------
    s : string
        String to attempt integer conversion.
    Return
    ------
    s : int or string
        Integer of input string if possible.  Otherwise string.

    """
    try:
        return int(s)
    except:
        return s


def natsort_key(s):
    """
    Used internally to get a list by which ``s`` is sorted.
    Called by :func:`natcmp`.

    Parameters
    -----------
    x :  string
        String to sort.
    Returns
    --------
    sort_key : list
        List used to sort input string.

    """
    sort_key = map(try_int, re.findall(r'(\d+|\D+)', s))
    return sort_key


def natcmp(a, b):
    """
    Natural string comparison. Case sensitive.

    Acts as replacement to basic Python standard library ``cmp`` for
    :func:`natsort`.

    Called by: :func:`natcasecmp`, :func:`natsort`,
    :func:`natsorted`.

    Parameters
    ----------
    a : string
        String to be compared with `b` (see Python stdlib `cmp` docs).
    b : string
        String to be compared with ``a``.

    Returns
    --------
    cmp_res : int
        Negative if ``a < b``, zero if ``a == b``, positive if ``a > b``.

    """
    return cmp(natsort_key(a), natsort_key(b))


def natcasecmp(a, b):
    """
    Natural string comparison. Ignores case.

    Acts as replacement to basic Python standard library ``cmp`` func for
    :func:`natsort`.

    Internal function.
    Can act as alternative to :func:`natcmp` with case insensitivity.

    Parameters
    ----------
    a : string
        String to be compared with ``b`` (see Python stdlib ``cmp`` docs).
    b : string
        String to be compared with ``a``.

    Returns
    --------
    cmp_res : int
        Negative if ``a < b``, zero if ``a == b``, positive if ``a > b``.

    """
    cmp_res = natcmp(a.lower(), b.lower())
    return cmp_res


def natsort(seq, cmp=natcmp):
    """
    Sorts input list or tuple *in place*.

    Parameters
    ----------
    seq : list or tuple
        List of strings to be sorted.
    cmp : function, [default = :func:`natcmp`]
        Function used to perform comparison.

    """
    seq.sort(cmp)


def natsorted(seq, cmp=natcmp):
    """
    Sorts input list or tuple as in :func:`natsort`, but returns a *copy*
    of original data rather than sorting in place.

    Parameters
    -----------
    seq : list or tuple
        List of strings to be sorted.
    cmp : function [default=:func:`natcmp`]
        Function used to perform comparison.

     Returns
     -------
     temp : list
        A copy of ``seq``, sorted by natural string sort :func:`natsort`.

    """
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp


#=============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    main()
