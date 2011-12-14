# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Provide basic utility functions for MD fluid simulation processing.

Most functions are of general use and could be used in other applications.

**Credits**: *Common string prefix*: Ned Batchelor, Cog project.

"""

# TODO: Ability to handle greater variety of filenames via regex.
# TODO: Add function which can take general user input subject to useer defined
#       constraints.  Could be basis for more specific existing yes/no question
#       in user_approve(...).


#-----------------------------------------------------------------------------#

import os
import re
import copy
import glob

from mdhandle.logger import Logger
from mdhandle.natsort  import natsort

#-----------------------------------------------------------------------------#

logger = Logger()

#-----------------------------------------------------------------------------#


def file_list(fn_base, verbose=False):
    """
    Generate a list of files matching glob of ``fn_base``.

    Parameters
    -----------
    fn_base : string
        Absolute path to target files.  This is input to glob expression.
    verbose : boolean
        Flag for verbosity of user input/output.

    Returns
    --------
    fglob : list
        List of files matching glob operation using ``fn_base``.

    """
    if fn_base is None:
        fn_base = ''

    fglob = glob.glob(fn_base)

    if (len(fglob) < 1) and (verbose is True):
        logger.warning('Did not find one or more files for basename given.')

    # sort flist inplace "naturally" in place: function defined below
    natsort(fglob)

    return fglob


def cull_flist_by_time(flist, time_list, verbose=False):
    """
    Eliminate entries within list of files which do not match ``time_list``.

    Parameters
    ----------
    flist : iterable
        List of filenames.
    time_list : iterable
        List of allowable times.

    Returns
    --------
    fn_return: list
        Sorted version of ``flist`` with times matching ``time_list``.
    fn_culled : list
        List of files with times falling outside of ``time_list``.

    """
    # TODO: replace with slice operator is possible and doesn't result in an
    # empty list being returned

    if hasattr(flist, '__iter__') is not True:
        if verbose is True:
            logger.warning('List of times must be a python list.')
        return []

    if len(flist) == 0 or len(time_list) == 0:
        return []

    # sort flist inplace "naturally" in place: function defined below
    natsort(flist)

    # time_list in ascending order
    time_list.sort()

    fn_return = []
    fn_culled = []
    time_remaining = copy.deepcopy(time_list)

    for fn in flist:
        time = time_from_filename(os.path.basename(fn))

        if time in time_list:
            fn_return.append(fn)
            time_there = True
            while time_there:
                try:
                    time_remaining.pop(time_remaining.index(time))
                except ValueError:
                    time_there = False
        else:
            fn_culled.append(fn)
            continue

    if verbose is True:
        if len(fn_culled) > 0:
            logger.user_message('The following list of files has been culled:')
            print_file_list(fn_culled)
        else:
            logger.user_message('No files removed from original list.')

        if len(time_remaining) > 0:
            logger.user_message('Times not found in the file list:')
            for t in time_remaining:
                logger.user_message(t)
        else:
            logger.user_message('All requested times were found in file list.')

    return fn_return, fn_culled


def cull_flist_by_function(flist, culling_function, verbose=False):
    """
    Culls files that do not exist from input file list.

    Parameters
    ----------
    flist : iterable
        List of filenames.
    culling_function : callable
        Callable which takes a filename as input and returns a boolean
        value {True | False}.  Must take a single input parameter.
    verbose:            boolean
        If ``True``, verbose user input and output.

    Returns
    -------
    fn_return : list
        Reduced list of files which return ``True`` when submitted to
        ``culling_function()``.
    fn_culled : list
        List of files for which ``culling_function`` returns ``False``.

    """
    if len(flist) == 0:
        return []

    fn_return = [f for f in flist if (culling_function(f) is True)]
    fn_culled = [f for f in flist if (culling_function(f) is False)]

    if verbose is True:
        if len(fn_culled) > 0:
            logger.user_message('The following list of files has been culled:')
            print_file_list(fn_culled)
        else:
            logger.user_message('No files removed from original list.')
    return fn_return, fn_culled


def time_from_filename(flist):
    """
    Extract the timestep value from a filename.

    Assume file name has a specific format with timestep in filename.
    (e.g. ``dump_100100.out``)

    Can also be used more generally to extract number from a string.

    Parameters
    -----------
    flist : iterable
        List of files with common basename and varying timesteps.

    Returns
    --------
    times : list
        List of integer timestep values.

    """
    # TODO: More flexible way to extract timestep from filename:
    #       - Default and user supplied regex (in settings.py?)

    # Check if flist is iterable
    if hasattr(flist, '__iter__'):
        times = []
        for fn in flist:
            # RECURSION in case have nested iterables
            times.append( time_from_filename(fn) )

    # Handles lowest level of recursion and case when string is given for flist
    elif isinstance(flist, basestring):
        fn = flist

        # HARDCODED: Currently assuming that numbers right before extension.
        times = int(re.findall(r'(\d+\.)', fn)[0][:-1])

    # RECURSION: For all other data types - convert to string and then recurse.
    else:
        fn = str(flist)
        times = time_from_filename(fn)

    return times


def print_file_list(flist, name=None):
    """
    Helper function for printing file lists.

    Parameters
    ------------
    flist : iterable
        List of file names.
    name : string, [default=None]
        Name of file list for header on printed output.

    """

    if name is not None:
        logger.user_message('Files contained in %s:' % name)
    # Do nothing
    else:
        pass

    for fn in flist:
        logger.user_message('\t %s' % fn)


def user_approve(question=''):
    """
    Asks user for yes/no approval.

    Parameters
    ----------
    question : string, [default='']
        Question for the user.

    Returns
    --------
    return_val : boolean
        Boolean ``True`` for 'yes' and boolean ``False`` for 'no'.

    """
    input_ok = False
    while not(input_ok):
        return_val = False
        
        response = raw_input(question + ' [Y/n]: ')
        response = response.lower()

        if response == 'y':
            return_val = True
        elif response == 'yes':
            return_val = True
        elif response == 'n':
            return_val = False
        elif response == 'no':
            return_val = False
        else:
            print '\nRespond with y or n to proceed.'
            continue
        return return_val


def user_python_input(question='', input_type=str):
    """
    Asks user for non yes/no question and returns user response, which is
    assumed to be a basic python data type.

    Parameters
    -----------
    question : string, [default='']
        Question for user.
    input_type : object, [default=str]
        Expected date type of user response.

    Returns
    ---------
    response : object
        Returns user response in the format of ``input_type``

    """
    input_ok = False
    while not(input_ok):
        response = input('%s type: %s' % (question, input_type))

        if isinstance(response, input_type) is True:
            input_ok = True
        else:
            input_ok = False
            print('User response does not match requsted data type')

    return response


def common_prefix(strings):
    """
    Finds the longest string that is a prefix of all the strings in supplied
    list.

    Credit: Ned Batchelor, Cog project.

    Parameters
    ----------
    strings : iterable
        List or tuple of strings.

    Returns
    -------
    prefix : string
        Common prefix shared by all strings in supplied list.

    """
    if not strings:
        return ''

    prefix = strings[0]
    for item in strings:
        if len(item) < len(prefix):
            prefix = prefix[:len(item)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != item[i]:
                prefix = prefix[:i]
                break
    return prefix


def clean(data, condition, fill_val=0.):
    """
    Replaces bad values within :class:`numpy.ndarray` with ``fill_val``
    based on user defined condition.

    Input ``data`` is modified in place.

    Useful for removing ``numpy.inf`` or ``numpy.NaN``.

    For two-dimensional arrays (vector and tensor collections),
    ``fill_val`` remains a simple object (i.e. single int or float)
    as replacement happens at element level.

    Parameters
    ------------
    data : class:`numpy.ndarray`
        :mod:`Numpy` data array.
        Shape and data type are arbitrary.
    condition : boolean
        Boolean array with same shape as data.
        Can be generated in place by boolean expression.
    fill_val : object
        Replaces entries in data if condition is ``True``.

    Returns
    --------
    data : :class:`numpy.ndarray`
        Cleaned version of ``data``.

    Examples
    ---------

    >>> # Replacing non-finite values with 0.
    >>> utilities.clean(a_vec, ~numpy.isfinite(s), 0.)

    """
    data[condition] = fill_val
    return data


#=============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    main()
