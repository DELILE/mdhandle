# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Provides functions for mapping column numbers, names, data types from LAMMPS
dump files.

Used in creating columnMapping metadata within HDF files.

"""

# TODO: Enable vector/tensors/symm tensors to read from file
# TODO: Abstraction for input testing to avoid repeated code.
# TODO: Merge number and no-number functions (argument)
# TODO: Merge user and file based results ( default args ? )

import sys
import os

import mdhandle.utilities as util
from mdhandle.logger import Logger
from mdhandle import settings

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


def lammps_mapping(dump_fn):
    """
    Reads column mapping of self-describing LAMMPS dump file from the 9th line.
    LAMMPS dumpe files were made self describing after approx 2007.

    Intended for use in reading LAMMPS ASCII files.

    Parameters
    ----------
    dump_fn : string
        Name of LAMMPS dump file name to be read.

    Returns
    --------
    columnsDict : dict
        Mapping of column number, name, type and location.

    """

    f = open(dump_fn, 'r')

    # Burning off all but the last line of the header.
    for i in range(len(settings.LAMMPSHEADER) - 1):
        f.readline()

    self_description = f.readline()

    # First two items are 'ITEM' and 'ATOMS'
    self_description = self_description.split()[2:]

    logger.user_message('Input dtype [int | float ] and')
    logger.user_message('location [Node | Cell]:')


    columnsDict = {}
    for col, name in enumerate(self_description):
        goodInput = False
        while not goodInput:
            usrIn = raw_input("%s %s --> " % (col, name))
            if len(usrIn.split()) != 2:
                logger.bad_input('Input format is: dtype location')
                continue
            else:
                dtype = str(usrIn.split()[0]).lower()
                location = str(usrIn.split()[1]).title()

            if location not in settings.LOCATIONS:
                logger.bad_input("Data location must be %s." %
                                                            settings.LOCATIONS)
                continue

            if dtype not in settings.DTYPES:
                logger.bad_input("dtype must be %s." % settings.DTYPES)
                continue

            columnsDict[name] = dict(col=col, dtype=dtype, location=location)
            goodInput = True

    logger.user_message('Scalar columns:')
    print_ordered_columns(columnsDict)

    if util.user_approve('Is this mapping correct?') is True:
        if util.user_approve('Print column data to file?') is True:
            out_fn = logger.request('Column file name: ')
            outf = open(out_fn, 'w')
            outf.write('# column name dtype location\n')
            for k in columnsDict:
                outf.write("%d %s %s %s\n" % (columnsDict[k]['col'], k,
                                              columnsDict[k]['dtype'],
                                              columnsDict[k]['location']))
        return columnsDict
    else:
        logger.bad_input('Going to user input for column mapping.')
        return user_mapping(self_description)


def file_mapping(col_file='', required=None, comment='#'):
    """
    Reads file mapping column number to names,  data type and grid location for
    reading LAMMPS ASCII data files.

    Fully describes the construction of all scalars from LAMMPS dump file.

    File format: number name dtype location

    Parameters
    ----------
    col_file : string, [default='']
        File to be used for mapping names and types to columns.
        Must be in format ``'n name type location'``.
        
        Type must be {'Int' | 'Float'}.
        
        Location is optiontal, but must be {'Node' | 'Cell'}, 
        [Default = 'Node']

    required : iterable, [default=None]
        List of required column names.
    comment : string, [default='#']
        Comment character

    Return
    -------
    columnsDict : dict
        Mapping of column names, dtypes, locations to column number.

    """
    columnsDict = {}

    if os.path.exists(col_file) is not True:
        logger.bad_input('Column file provided does not exist.')
        logger.bad_input('Going to user input')
        return user_mapping(required)

    f = open(col_file)

    for line in f:
        # Removing white space before and after.
        line = line.strip()

        # Ignoring header with '#' comment char.
        if line[0] == comment:
            continue

        n = int(line.split()[0])
        name = str(line.split()[1])
        dtype = str(line.split()[2]).lower()
        location = str(line.split()[3]).title()

        if location not in settings.LOCATIONS:
            logger.bad_input("Data location must be %s." % settings.LOCATIONS)
            user_mapping(required)

        if dtype not in settings.DTYPES:
            logger.bad_input("Data type must be %s." % settings.DTYPES)
            user_mapping(required)

        columnsDict[name] = dict(col=n, dtype=dtype, location=location)

    f.close()
    
    if required is None:
        pass
    elif all( [(i in columnsDict) for i in required] ) is not True:
        logger.bad_input('%s required. Restarting' % required)
        user_mapping(required)

    logger.user_message('Scalar columns:')
    print_ordered_columns(columnsDict)

    if util.user_approve('Is this mapping correct?') is True:
        return columnsDict
    else:
        logger.bad_input('Restarting user input for column mapping.')
        return user_mapping(required)


def user_mapping(required=None):
    """
    Takes user input to identify the column name, number, data type
    and grid location in LAMMPS ASCII data files.

    * User input must be in format ``'n name type location'``.
    * No spaces alled within name
    * Type must be {'int' | 'float'}
    * Grid location must be {'Node' | 'Cell'}

    Parameters
    -----------
    required : iterable, [default=None]
        List of required data columns.

    Returns
    -------
    columnsDict : dict
       Mapping of column names, dtypes, locations to column number.

    """
    logger.procedure_banner('User Input for Scalar Identification')

    logger.request("Input number, name, dtype and location\n\
                    (zero-indexed, blank line to end): ")

    usrIn = None
    columnsDict = {}

    while usrIn != '':
        usrIn = raw_input("--> ")

        # loop exit condition
        if usrIn == '':
            continue

        if len(usrIn.split()) != 4:
            logger.bad_input('Input format is: number name dtype location')
            continue

        try:
            n = int(usrIn.split()[0])
        except ValueError:
            logger.bad_input("Column number must be an integer")
            continue
        name = str(usrIn.split()[1])
        dtype = str(usrIn.split()[2]).lower()
        location = str(usrIn.split()[3]).title()

        if location not in settings.LOCATIONS:
            logger.bad_input("Data location must be %s." % settings.LOCATIONS)
            continue

        if dtype not in settings.DTYPES:
            logger.bad_input("dtype must be %s." % settings.DTYPES)
            continue

        columnsDict[name] = dict(col=n, dtype=dtype, location=location)

    if required is None:
        pass
    elif all( [ (i in columnsDict) for i in required]) is not True:
        logger.bad_input('%s required. Restarting' % required)
        user_mapping(required)

    logger.user_message('Scalar columns:')
    print_ordered_columns(columnsDict)

    if util.user_approve('Is this mapping correct?') is True:
        if util.user_approve('Print column data to file?') is True:
            logger.request('Column file name: ')
            out_fn = raw_input('>>> ')
            outf = open(out_fn, 'w')
            for k in columnsDict:
                outf.write("%d %s %s %s\n" % (columnsDict[k]['col'], k,
                                              columnsDict[k]['dtype'],
                                              columnsDict[k]['location']))
        return columnsDict
    else:
        logger.bad_input('Restarting user input for column mapping.')
        return user_mapping(required)


def convert_from_num(col_mapping):
    """
    Removes numbering in column mapping.

    Parameters
    ----------
    col_mapping : dict
        Column mapping linking item name to column number, data type and
        location.

    Returns
    --------
    columnsDict : dict
        Returns column mapping without column number.

    """
    return dict( (k, dict(dtype=col_mapping[k]['dtype'],
                          location=col_mapping[k]['location'] )  )
                   for k in col_mapping  )


def _identify_collections(col_mapping, kind='vectors'):
    """
    Produces dictionary containing the name (key),  components, data type
    and grid location (i.e. 'Node', 'Cell') of vectors,tensors and symmetric
    tensors.

    Based on user input at the command line.
    Intended for use in one time conversion of LAMMPS output.

    Parameters
    ----------
    col_mapping : dict
        Dictionary mapping scalar name, column number, dtype and location.
    kind : string, {'vectors' | 'tensors' | 'symm_tensors'}, [default='vectors']
        Type of data to be identified.

    Returns
    --------
    collection : dict
        Returns dictionary mapping name of collection to components, data
        type and grid location.

    """
    if kind not in settings.PROP_SIZE and kind != 'scalars':
        logger.error('Invalid collection type.  \n\
                      \t Must be: "vector", "tensor", or "symm_tensor".')

    # ----------------- USER INPUT  -------------------------------------------

    logger.procedure_banner('User Input for %s Identification' % kind)
    logger.user_message('Input: name and components')
    logger.user_message('(space separated, blank line to end):')

    usrIn = None
    collection = {}

    while usrIn != '':
        usrIn = logger.request('', input_type='raw')

        # loop exit condition
        if usrIn == '':
            continue

        elif len(usrIn.split()) != (1 + settings.PROP_SIZE[kind]):
            logger.bad_input("Input must be at least %d words."
                                            % (1 + settings.PROP_SIZE[kind]))
            continue

        # Good input
        else:
            name = str(usrIn.split()[0])
            components = [str(i) for i in usrIn.split()[1:]]

        # Looking for spaces in name and components
        spaces = ' ' in name  or  any([' ' in i for i in components])

        if spaces is True:
            logger.bad_input('No spaces in attribute name.')
            continue

        if all([(i in col_mapping) for i in components]) is not True:
            logger.bad_input('Components must exist in col_mapping.')
            continue

        if all([col_mapping[i]['dtype'] == col_mapping[components[0]]['dtype']]
                                        for i in components ):
            dtype = col_mapping[components[0]]['dtype']
        else:
            logger.bad_input('dtype must be homogeneous.')
            continue

        if all([col_mapping[i]['location'] ==
                  col_mapping[components[0]]['location']] for i in components):
            location = col_mapping[components[0]]['location']
        else:
            logger.bad_input('Location of components must be homogeneous.')
            continue

        collection[name] = dict(comps=components, dtype=dtype,
                                                  location=location)

    if collection == {}:
        logger.user_message('No %s are specified.' % kind)
    else:
        logger.user_message('%s are: ' % kind.title())
        for k in collection:
            logger.user_message('%s --> %s %s %s' %
                                    (k, collection[k]['comps'],
                                        collection[k]['dtype'],
                                        collection[k]['location'] ) )

    if util.user_approve('Is this list of %s correct?' % kind ) is True:
        logger.completion_msg('COMPLETED %s IDENTIFICATION' % kind)
        return collection
    else:
        logger.user_message('Returning to user input for %s.' % kind)
        return _identify_collections(col_mapping, kind)


def vector_mapping(col_mapping):
    """
    Identifies vector collections within existing scalar data.

    Wrapper to :func:`_identify_collections`.

    """
    return _identify_collections(col_mapping, 'vectors')


def tensor_mapping(col_mapping):
    """
    Identifies tensor collections within existing scalar data.

    Wrapper to :func:`_identify_collections`.

    """
    return _identify_collections(col_mapping, 'tensors')


def symm_tensor_mapping(col_mapping):
    """
    Identifies symmetric tensor collections within existing scalar data.

    Wrapper to :func:`_identify_collections`.

    """
    return _identify_collections(col_mapping, 'symm_tensors')


def print_ordered_columns(col_mapping):
    """
    Take a numbered ``col_mapping`` dictionary and print in ordered way,
    based on column number.

    Will fail if entries in ``col_mapping`` dict does not contain ``'col'``.

    Parameters
    ----------
    col_mapping : dict
        Dictionary with column number, dtype and location stored under name.

    """
    # Sorting by column numbers.
    sorted_cols = [(k, v) for (k, v) in col_mapping.items()]
    try:
        sorted_cols.sort(key=lambda item: item[1]['col'])
    except KeyError:
        # Catching case where col_mapping items do not have 'col'
        logger.warning('Column mapping does not have column number.')
        return

    for item in sorted_cols:
        logger.user_message('%s --> %s %s %s' % (item[0],
                                                 item[1]['col'],
                                                 item[1]['dtype'],
                                                 item[1]['location'] ) )


#==============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
