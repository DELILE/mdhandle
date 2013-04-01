# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Manages XDMF files used to contain snapshot data from MD simulation including
atomistic and gridded data..

Produces required ``.xmf`` and HDF files for a single snapshot.

Updates existing ``.xmf`` and HDF files for single snapshot

Creates summary ``.xmf`` files for temporal collection of snapshots

"""
#TODO: xml parser/writer to handle xmf action.
#TODO: Use XDMF library: Add parallel IO functionality from XDMF library

#TODO: Add function at module level to ouput a spatial summary file, which
#      brings together disparate spatial grids (as opposed to separate xmf
#      files for different timesteps).
# TODO: Add HDF function to move array from one data set to another
# TODO: change funciton/method names so more abstract (i.e. another writer
#       could be used)
# TODO: Wrapper around file opening to keep count of the number of open files.
# TODO: Migrate xmf scalar, vector and tensor to just wrap xmf_block_attribute

# TODO: Add function to check if file is open before write, OR handle file is closed exceptions gracefully.

#----------------------------------------------------------------------------

import os
import shutil

import tables as pyT

from mdhandle import utilities as util
from mdhandle.logger import Logger
from mdhandle import vectors
from mdhandle import settings
from mdhandle.properties.binned_properties import IMPLEMENTED_GRIDS

#----------------------------------------------------------------------------
#   Module functions
#----------------------------------------------------------------------------

logger = Logger()

#----------------------------------------------------------------------------


def xmf_write_temporal_summary(sim_cont, fn_out='TEMP-SUMMARY-DEFAULT.xmf'):
    """
    Creates a summary xmf file from set of individual snapshot xmf files.
    Assumes each snapshot xmf file contains a single uniform grid.
    The file may be for raw simulation data or may be the result of
    processing.

    XML content within and including grid is copied into summary xmf file.
    
    .. note::
       XMF temporal summary file must be in the same directory as the 
       individual XMF files.

    Parameters
    ----------
    sim_cont : :class:`mdhandle.simcontainer.SimContainer`
    fn_out : string,  [default=SUMMARY-DEFAULT.xmf]
        FileName for summary file.  Do not give absolute path.

    Returns
    -------
    0
        Successful execution.
    -1
        An error has occured and the XMF file may not have been generated
        at all.

    """
    if sim_cont.meta['grid_type'] not in IMPLEMENTED_GRIDS:
        logger.error('Grid type %s unknown' % sim_cont.meta['grid_type'])

    if fn_out == 'TEMP-SUMMARY-DEFAULT.xmf':
        # Handling the grid_type - changes the resulting file name
        if sim_cont.meta['grid_type'] == 'atoms':
            add_name = 'TEMPORAL-SUMMARY'
        elif sim_cont.meta['grid_type'] == '3DRegular':
            add_name = sim_cont.active_dataset_name + '_TEMPORAL-SUMMARY'

        fn_out = '%s%s.xmf' % (util.common_prefix(sim_cont.flist), add_name)

    text = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.1">
    <Domain>
        <Grid Name="XXX_TEMPGRIDNAME_XXX" GridType="Collection" CollectionType="Temporal">
"""
    temporal_grid_name = '%s-%s-%s-Temporal-Collection' %\
                            (sim_cont.meta['sim_name'],
                             sim_cont.meta['grid_type'],
                             sim_cont.active_dataset_name)

    text = text.replace("XXX_TEMPGRIDNAME_XXX", temporal_grid_name)

    for snap in sim_cont.slist:
        xmf_fn = snap.writer.xmf_snap_name()
        text += '           <xi:include href="%s" xpointer="xpointer(//Xdmf/Domain/Grid)" />\n' % os.path.basename(xmf_fn)

    # Adding final text and writing to file
    text += """
    </Grid>
</Domain>
</Xdmf>"""

    if os.path.exists(fn_out) is True:
        logger.warning('Overwriting existing XMF summary file: %s' % fn_out)
        os.remove(fn_out)
    
    f_out = open(fn_out, 'w')
    f_out.write(text)
    f_out.close()
    return 0


def hdf_compact(flist, verbose=False):
    """
    Performs ptrepack compaction on list of HDF files.  This procedure
    needs to be done periodically after moving data around within the HDF
    tree to prevent the file size from increasing.

    flist : iterable
        List of HDF files to be compacted.
    verbose : boolean,  [default=False]
        Flag for verbosity of user input and output.

    Returns
    -------
    0
        Successfully completed compaction.
    -1
        An error has occured.  Unable to compact some files.

    """
    # ASSERT flist is iterable
    assert hasattr(flist, '__iter__')

    error = False
    error_files = []

    # TODO: hard coded tmp directory for all files - change to setting.
    for fn in flist:
        fn_tmp = os.path.join('/tmp', fn+'_tmp.h5')
        success = os.system('ptrepack -o %s %s' % (fn, fn_tmp))
        if success == 0:
            shutil.move(fn_tmp, fn)
        else:
            error = True
            error_files.append(fn)
            os.remove(fn_tmp)

    if error is True:
        if verbose is True:
            logger.warning('Compaction Failed for: ')
            for errorFile in error_files:
                logger.warning('    %s' % os.path.basename(errorFile))
        return -1
    else:
        logger.completion_msg('Successfully compacted set of HDF files.')
        return 0


def valid_file(fn):
    """
    Tests whether ``fn`` is a valid HDF file.

    Parameters
    -----------
    fn : string
        Pathname to snapshot file.

    Returns
    -------
    is_valid : boolean
        ``True`` if ``fn`` is a valid HDF file that can be read by the
        Py-Tables package.

    """
    # Can add other conditions as needed through logical and
    is_valid = pyT.isHDF5File(fn)
    return is_valid


#-----------------------------------------------------------------------------
# Classes
#-------------------------------------------------------------------------------


class XDMFWriter(object):
    """
    Creates and/or maintains existing HDF files for storing snapshot data.

    Creates new ``.xmf`` files.

    Parameters
    -----------
    snap :   :class:`mdhandle.snap.Snap`
        Simulation snapshot.
    verbose : boolean,  [default=False]
        Flag for verbosity of user input and output.

    Attributes
    -----------
    snap : :class:`mdhandle.snap.Snap`
    logger : :class:`mdhandle.logger.Logger`
        Logger and user I/O module.
    f_hdf : ``HDF file``
        HDF file handle opened by PyTables.
    verbose : boolean
        If ``True``, verbose user I/O.
    _xmf : file
        Handle for ``.xmf`` file.

    """

    def __init__(self, snap, verbose=False):
        self.snap = snap
        self.verbose = verbose
        self.logger = Logger()
        self.f_hdf = None
        self._xmf = ''


    # -------------------------------------------------------------------------
    #   XMF Methods.
    # -------------------------------------------------------------------------

    def xmf_snap_name(self):
        """
        Generates name of XMF file for :class:`mdhandle.snap.Snap`
        object.

        Returns
        --------
        xmf_name : string
            Properly formatted ``.xmf`` file name for single snap shot
            based on snap file name and dataset name.

        """
        if self.snap.meta['grid_type'] not in IMPLEMENTED_GRIDS:
            logger.error('Grid type %s unknown' % self.snap.meta['grid_type'])

        if self.snap.meta['grid_type'] == 'atoms':
            xmf_name = os.path.splitext(self.snap.filename)[0]+'.xmf'

        # TODO: Recast in terms fo snap attributes (e.g. sim_name, etc.)
        elif self.snap.meta['grid_type'] == '3DRegular':
            # TODO: HARDCODED - Assume format 'dump_simname_1234.h5'
            splits = os.path.splitext(self.snap.fn_base)[0].split('_')
            leader = '_'.join(splits[:-1])
            xmf_name = '%s_%s_%s.xmf' % (leader, self.snap.active_dataset_name,
                                                 self.snap.meta['time'])
            xmf_name = os.path.join(self.snap.directory, xmf_name)
        return xmf_name

    def xmf_new(self):
        """
        Creates a new xmf file.

        """
        self._xmf = open(self.xmf_snap_name(), 'w')

        XML = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.1">
    <Domain>"""
        self._xmf.write(XML)

    def xmf_finish(self):
        """
        Writes final XML to ``.xmf`` file, and closes ``.xmf`` file.

        """

        XML = """
    </Domain>
</Xdmf>"""
        self._xmf.write(XML)
        self._xmf.close()

        if self.verbose:
            self.logger.completion_msg('Wrote metadata to %s' % 
                                                self.xmf_snap_name())

    def xmf_grid_init(self):
        """
        Writes initial generic XML for both atom and gridded data
        to ``.xmf`` file.

        """
        grid_name = '%s-%s-%s' % (self.snap.meta['sim_name'],
                                  self.snap.meta['grid_type'],
                                  self.snap.active_dataset_name)
        XML = """
            <Grid Time="XXX_TIME_XXX" Name="XXX_NAME_XXX" GridType="Uniform">
                <Time Value="XXX_TIME_XXX" />"""
        XML = XML.replace("XXX_TIME_XXX", str(self.snap.meta['time']))
        XML = XML.replace("XXX_NAME_XXX", str(grid_name))
        self._xmf.write(XML)

    def xmf_grid_atoms(self, num_atoms):
        """
        Writes XML for topology and geometry information for unstructured
        ``Polyvertex`` grid.

        .. note::
           Does not write ``xyz`` position of grid points to XML.
           These are stored in the HDF file.

        Parameters
        -----------
        num_atoms : int
            Number of atoms in grid.

        """

        XML = """
            <Topology TopologyType="Polyvertex" NumberOfElements="XXX_NELEMENTS_XXX"/>
            <Geometry GeometryType="X_Y_Z">
                <DataItem Dimensions="XXX_NPOINTS_XXX" NumberType="Float" Format="HDF">
                    XXX_FBASE_XXX:/XXX_DATASET_XXX/x
                </DataItem>
                <DataItem Dimensions="XXX_NPOINTS_XXX" NumberType="Float" Format="HDF">
                    XXX_FBASE_XXX:/XXX_DATASET_XXX/y
                </DataItem>
                <DataItem Dimensions="XXX_NPOINTS_XXX" NumberType="Float" Format="HDF">
                    XXX_FBASE_XXX:/XXX_DATASET_XXX/z
                </DataItem>
            </Geometry>"""

        XML = XML.replace("XXX_NELEMENTS_XXX", str(num_atoms))
        XML = XML.replace("XXX_NPOINTS_XXX", str(num_atoms))
        XML = XML.replace("XXX_FBASE_XXX", self.snap.fn_base)
        XML = XML.replace("XXX_DATASET_XXX", self.snap.active_dataset_name)
        self._xmf.write(XML)

    def xmf_grid_regular(self, spacing, origin, num_elem):
        """
        Writes XML for topology and geometry for ``3DCORECTMesh``
        (i.e. regular 3D grid orthogonal to coordinate axes).

        Parameters
        -----------
        spacing : array_like
            Vector (``length == 3``) for grid spacing along ``[x,y,z]`` axes.
            Each element is either ``float`` or ``int``.
        origin : array_like
            Vector (``length == 3``) for location of origin of grid.
            Each element is either ``float`` or ``int``.
        num_elem :   array_like
            Vector (``length == 3``) for number of grid elements along
            ``[x,y,z]`` axes. Each element is ``int``.

        """

        XML = """
            <Topology TopologyType="3DCORECTMesh" NumberOfElements="XXX_NX_XXX XXX_NY_XXX XXX_NZ_XXX "/>
            <Geometry GeometryType="ORIGIN_DXDYDZ">
                <DataItem Dimensions="3 " NumberType="Float" Format="XML">
                    XXX_ORGX_XXX XXX_ORGY_XXX XXX_ORGZ_XXX
                </DataItem>
                <DataItem Dimensions="3 " NumberType="Float" Format="XML">
                    XXX_DX_XXX XXX_DY_XXX XXX_DZ_XXX
                </DataItem>
            </Geometry>"""
        num_elem = num_elem.astype('int')

        XML = XML.replace("XXX_NX_XXX", str(num_elem[0]))
        XML = XML.replace("XXX_NY_XXX", str(num_elem[1]))
        XML = XML.replace("XXX_NZ_XXX", str(num_elem[2]))

        XML = XML.replace("XXX_ORGX_XXX", str(origin[0]))
        XML = XML.replace("XXX_ORGY_XXX", str(origin[1]))
        XML = XML.replace("XXX_ORGZ_XXX", str(origin[2]))

        XML = XML.replace("XXX_DX_XXX", str(spacing[0]))
        XML = XML.replace("XXX_DY_XXX", str(spacing[1]))
        XML = XML.replace("XXX_DZ_XXX", str(spacing[2]))

        self._xmf.write(XML)

    def xmf_grid_finish(self):
        """
        Writes out finising XML Grid tag and closes ``.xmf`` file.

        """

        XML = """
            </Grid>"""
        self._xmf.write(XML)

    def xmf_scalar_attribute(self, name, dtype, location, length):
        """
        Writes XML for scalar attribute.
        Should only be used once grid/topology/geometry are already setup.

        Parameters
        ----------
        name : string
            Name of attribute in HDF file.
        dtype : string, {'Int' | 'Float'}
            Data type
        location : string,  {'Node' | 'Cell'}
            Location of scalar attribute within grid cell.
        length : int
            Number of scalars (e.g. number of atoms or grid points).

        """
        XML = """
                    <Attribute Type="Scalar" Center="XXX_ATTRLOC_XXX" Name="XXX_ATTRNAME_XXX">
                        <DataItem DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX" Format="HDF">
                            XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_ATTRNAME_XXX
                        </DataItem>
                    </Attribute>"""
        
        dtype = dtype.lower()
        if 'int' in dtype:
            dtype = 'Int'
        elif 'float' in dtype:
            dtype = 'Float'
        else:
            raise UserWarning('Invalid data type submitted. \
                                                    Must be int* or float*.')

        # Replace the markers
        XML = XML.replace("XXX_ATTRLOC_XXX", str(location))
        XML = XML.replace("XXX_ATTRNAME_XXX", str(name))
        XML = XML.replace("XXX_ATTRLEN_XXX", str(length))
        XML = XML.replace("XXX_DTYPE_XXX", str(dtype))
        XML = XML.replace("XXX_FBASE_XXX", self.snap.fn_base)
        XML = XML.replace("XXX_DATASET_XXX", self.snap.active_dataset_name)

        self._xmf.write(XML)

    def xmf_vector_attribute(self, name, dtype, location,
                                             length, v0_name, v1_name, v2_name):
        """
        Writes XML for vector attribute.
        Should only be used once grid/topology/geometry are already setup.

        Parameters
        -----------
        name : string
            Name of attribute in HDF file.
        dtype : string,  {'Int' | 'Float'}
            Data type.
        length : int
            Number of individual vectors within attribute
            (e.g  number of atoms or grid points).
        location : string,  {'Node' | 'Cell'}
            Location of scalar attribute within grid cell.
        v{0,1,2}_name : string
            Names of scalar values are to be the vector's components.

        """

        XML = """
                    <Attribute Type="Vector" Center="XXX_ATTRLOC_XXX" Name="XXX_ATTRNAME_XXX">
                        <DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="XXX_ATTRLEN_XXX 3">
                            <DataItem Name="XXX_V0_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_V0_XXX
                            </DataItem>
                            <DataItem Name="XXX_V1_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_V1_XXX
                            </DataItem>
                            <DataItem Name="XXX_V2_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_V2_XXX
                            </DataItem>
                        </DataItem>
                    </Attribute>"""
        dtype = dtype.lower()
        if 'int' in dtype:
            dtype = 'Int'
        elif 'float' in dtype:
            dtype = 'Float'
        else:
            self.logger.error('Invalid data type submitted. \
                                                      Must be int* or float*.')

        # Replace the markers
        XML = XML.replace("XXX_ATTRLOC_XXX", str(location))
        XML = XML.replace("XXX_ATTRNAME_XXX", str(name))
        XML = XML.replace("XXX_ATTRLEN_XXX", str(length))
        XML = XML.replace("XXX_DTYPE_XXX", str(dtype))
        XML = XML.replace("XXX_FBASE_XXX", self.snap.fn_base)
        XML = XML.replace("XXX_DATASET_XXX", self.snap.active_dataset_name)
        XML = XML.replace("XXX_V0_XXX", str(v0_name))
        XML = XML.replace("XXX_V1_XXX", str(v1_name))
        XML = XML.replace("XXX_V2_XXX", str(v2_name))

        self._xmf.write(XML)

    def xmf_symm_tensor_attribute(self, name, dtype, location, length,
                                    t0_name, t1_name, t2_name,
                                    t3_name, t4_name, t5_name):
        """
        Writes XML for symmetric tensor attribute.
        Should be used once grid/topology/geometry are already setup.

        Parameters
        -----------
        name : string
            Name of symmetric tensor attribute in HDF file.
        dtype : string,  {'Int' | 'Float'}
            Data type.
        location : string,  {'Node'|'Cell'}
            Location of scalar attribute within grid cell.
        length : int
            Number of individual tensors within attribute
            (e.g. number of atoms or grid points).
        t{0,1,2,3,4,5}_name : string
            Names of scalar values which are vector's components.
            Order is: ``00, 01, 02, 11, 12, 22``, which is the same as the
            order in LAMMPS (i.e. ``xx yy zz xy xz yz == 00 11 22 01 02 12``)

        """

        XML = """
                    <Attribute Type="Tensor6" Center="XXX_ATTRLOC_XXX" Name="XXX_ATTRNAME_XXX">
                        <DataItem ItemType="Function" Function="JOIN($0, $1, $2, $3, $4, $5)" Dimensions="XXX_ATTRLEN_XXX 6">
                            <DataItem Name="XXX_T0_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T0_XXX
                            </DataItem>
                            <DataItem Name="XXX_T1_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T1_XXX
                            </DataItem>
                            <DataItem Name="XXX_T2_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T2_XXX
                            </DataItem>
                            <DataItem Name="XXX_T3_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T3_XXX
                            </DataItem>
                            <DataItem Name="XXX_T4_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T4_XXX
                            </DataItem>
                            <DataItem Name="XXX_T5_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T5_XXX
                            </DataItem>
                        </DataItem>
                    </Attribute>"""
        dtype = dtype.lower()
        if 'int' in dtype:
            dtype = 'Int'
        elif 'float' in dtype:
            dtype = 'Float'
        else:
            self.logger.error('Invalid data type submitted. \
                                                      Must be int* or float*.')

        # Replace the appropriate markers
        XML = XML.replace("XXX_ATTRLOC_XXX", str(location))
        XML = XML.replace("XXX_ATTRNAME_XXX", str(name))
        XML = XML.replace("XXX_ATTRLEN_XXX", str(length))
        XML = XML.replace("XXX_DTYPE_XXX", str(dtype))
        XML = XML.replace("XXX_FBASE_XXX", self.snap.fn_base)
        XML = XML.replace("XXX_DATASET_XXX", self.snap.active_dataset_name)
        XML = XML.replace("XXX_T0_XXX", str(t0_name))
        XML = XML.replace("XXX_T1_XXX", str(t1_name))
        XML = XML.replace("XXX_T2_XXX", str(t2_name))
        XML = XML.replace("XXX_T3_XXX", str(t3_name))
        XML = XML.replace("XXX_T4_XXX", str(t4_name))
        XML = XML.replace("XXX_T5_XXX", str(t5_name))
        self._xmf.write(XML)

    def xmf_tensor_attribute(self, name, dtype, location, length,
                                t0_name, t1_name, t2_name,
                                t3_name, t4_name, t5_name,
                                t6_name, t7_name, t8_name):
        """
        Writes XML for general tensor attribute.
        Should be used once grid/topology/geometry are already setup.

        Parameters
        -----------
        name : string
            Name of tensor attribute in HDF file.
        dtype : string,  {'Int' | 'Float'}
            Data type.
        location : string,  {'Node' | 'Cell'}
            Location of scalar attribute within grid cell.
        length : int
            Number of individual tensors within attribute
            (e.g. number of atoms).
        t{0:8}_name: string
            Names of scalar values which are the tensor components. Order is:            
            ``00 01 02 10 11 12 20 21 22 or xx xy xz yx yy yz zx zy zz``

        """
        XML = """
                    <Attribute Type="Tensor" Center="XXX_ATTRLOC_XXX" Name="XXX_ATTRNAME_XXX">
                        <DataItem ItemType="Function" Function="JOIN($0, $1, $2, $3, $4, $5, $6, $7, $8)" Dimensions="XXX_ATTRLEN_XXX 9">
                            <DataItem Name="XXX_T0_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T0_XXX
                            </DataItem>
                            <DataItem Name="XXX_T1_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T1_XXX
                            </DataItem>
                            <DataItem Name="XXX_T2_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T2_XXX
                            </DataItem>
                            <DataItem Name="XXX_T3_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T3_XXX
                            </DataItem>
                            <DataItem Name="XXX_T4_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T4_XXX
                            </DataItem>
                            <DataItem Name="XXX_T5_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T5_XXX
                            </DataItem>
                            <DataItem Name="XXX_T6_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T6_XXX
                            </DataItem>
                            <DataItem Name="XXX_T7_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T7_XXX
                            </DataItem>
                            <DataItem Name="XXX_T8_XXX" Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX">
                                XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_T8_XXX
                            </DataItem>
                        </DataItem>
                    </Attribute>"""
        dtype = dtype.lower()
        if 'int' in dtype:
            dtype = 'Int'
        elif 'float' in dtype:
            dtype = 'Float'
        else:
            self.logger.error('Invalid data type submitted. \
                                                      Must be int* or float*.')

        # Replace the appropriate markers
        XML = XML.replace("XXX_ATTRLOC_XXX", str(location))
        XML = XML.replace("XXX_ATTRNAME_XXX", str(name))
        XML = XML.replace("XXX_ATTRLEN_XXX", str(length))
        XML = XML.replace("XXX_DTYPE_XXX", str(dtype))
        XML = XML.replace("XXX_FBASE_XXX", self.snap.fn_base)
        XML = XML.replace("XXX_DATASET_XXX", self.snap.active_dataset_name)
        XML = XML.replace("XXX_T0_XXX", str(t0_name))
        XML = XML.replace("XXX_T1_XXX", str(t1_name))
        XML = XML.replace("XXX_T2_XXX", str(t2_name))
        XML = XML.replace("XXX_T3_XXX", str(t3_name))
        XML = XML.replace("XXX_T4_XXX", str(t4_name))
        XML = XML.replace("XXX_T5_XXX", str(t5_name))
        XML = XML.replace("XXX_T6_XXX", str(t6_name))
        XML = XML.replace("XXX_T7_XXX", str(t7_name))
        XML = XML.replace("XXX_T8_XXX", str(t8_name))

        self._xmf.write(XML)

    def xmf_block_attribute(self, name, attribute_type,
                                    dtype, location, length, data_name):
        """
        Writes XML for data (scalars, vectors, tensors) stored as a
        contiguous block of data rather than in separate blocks for each
        component.

        Should be used once grid/topology/geometry are already setup.

        Parameters
        -----------
        name : string
            Name of attribute.
        attribute_type : string,  
            {'scalar' | 'vector' | 'symm_tensor' | 'tensor' }
            Type of attribute.
        dtype : string,  {'int' | 'float'}
            Data type.
        location : string,  {'Node' | 'Cell'}
            Location of scalar attribute within grid cell.
        attrLen : int
            Number of individual items within attribute
            (e.g. number of atoms).
        data_name : string
            Name of data in HDF file dataset.

        """

        XML = """
                    <Attribute Type="XXX_ATTRTYPE_XXX" Center="XXX_ATTRLOC_XXX" Name="XXX_ATTRNAME_XXX" Dimensions="XXX_ATTRLEN_XXX XXX_ATTRWIDTH_XXX">
                        <DataItem Format="HDF" DataType="XXX_DTYPE_XXX" Dimensions="XXX_ATTRLEN_XXX XXX_ATTRWIDTH_XXX">
                            XXX_FBASE_XXX:/XXX_DATASET_XXX/XXX_DATA_XXX
                        </DataItem>
                    </Attribute>"""

        # Identifying type of attribute
        attribute_type = attribute_type.lower()
        if attribute_type == 'scalar':
            width = 1
        elif attribute_type == 'vector':
            width = 3
        elif attribute_type == 'symm_tensor':
            attribute_type = 'Tensor6'
            width = 6
        elif attribute_type == 'tensor':
            width = 9
        else:
            self.logger.error('Invalid XDMF type submitted must be scalar, \
                                                vector, symm_tensor or tensor')
        attribute_type = attribute_type.title()

        dtype = dtype.lower()
        if 'int' in dtype:
            dtype = 'int'
        elif 'float' in dtype:
            dtype = 'float'
        else:
            self.logger.error('Invalid data type submitted. \
                                                      Must be int* or float*.')

        # Replace the appropriate markers
        XML = XML.replace("XXX_ATTRTYPE_XXX", str(attribute_type))
        XML = XML.replace("XXX_ATTRLOC_XXX", str(location))
        XML = XML.replace("XXX_ATTRNAME_XXX", str(name))
        XML = XML.replace("XXX_ATTRLEN_XXX", str(length))
        XML = XML.replace("XXX_ATTRWIDTH_XXX", str(width))
        XML = XML.replace("XXX_DTYPE_XXX", str(dtype))
        XML = XML.replace("XXX_FBASE_XXX", self.snap.fn_base)
        XML = XML.replace("XXX_DATASET_XXX", self.snap.active_dataset_name)
        XML = XML.replace("XXX_DATA_XXX", str(data_name))

        self._xmf.write(XML)

    def xmf_snap_file(self):
        """
        Create a new xmf file.  This is light operation as works with
        metadata.

        """

        self.xmf_new()
        self.xmf_grid_init()

        if self.snap.meta['grid_type'] == 'atoms':
            self.xmf_grid_atoms(self.snap.meta['num_atoms'])
            data_len = self.snap.meta['num_atoms']
        elif self.snap.meta['grid_type'] == '3DRegular':
            self.xmf_grid_regular(self.snap.meta['grid_spacing'],
                                  self.snap.sim_cell.origin,
                                  self.snap.meta['num_elem'])
            data_len = ' '.join([str(i) for i in self.snap.meta['num_elem']])

        for (key, val) in self.snap.col_mapping.iteritems():
            if key in ( self.snap.vectors.keys() +
                        self.snap.tensors.keys() +
                        self.snap.symm_tensors.keys() ):
                continue
            self.xmf_scalar_attribute(key, val['dtype'], val['location'],
                                                         data_len)

        for (key, val) in self.snap.vectors.iteritems():
            if 'comps' in val:
                self.xmf_vector_attribute(key, val['dtype'], val['location'],
                                               data_len, *val['comps'])
            # Contiguous vector
            else:
                self.xmf_block_attribute(key, 'vector', val['dtype'],
                                              val['location'], data_len, key)

        for (key, val) in self.snap.tensors.iteritems():
            if 'comps' in val:
                self.xmf_tensor_attribute(key, val['dtype'], val['location'],
                                               data_len, *val['comps'])
            # Contiguous tensor
            else:
                self.xmf_block_attribute(key, 'tensor', val['dtype'],
                                            val['location'], data_len, key)

        for (key, val) in self.snap.symm_tensors.iteritems():
            if 'comps' in val:
                self.xmf_symm_tensor_attribute(key, val['dtype'],
                                val['location'], data_len, *val['comps'])
            # Contiguous symm_tensor
            else:
                self.xmf_block_attribute(key, 'symm_tensor', val['dtype'],
                                            val['location'], data_len, key)

        self.xmf_grid_finish()
        self.xmf_finish()

        if self.verbose:
            self.logger.completion_msg('Done making new/updated xmf for \
               %s dataset in %s' % (self.snap.active_dataset_name, 
                                    self.snap.fn_base) )

    # -------------------------------------------------------------------------
    #       HDF Functions
    # -------------------------------------------------------------------------

    def is_open(self):
        """
        If HDF file is open, returns ``True``.
        
        If HDF has not yet been opened at all, returns False.

        """
        if self.f_hdf is None:
            return False
        if self.f_hdf.isopen == 1:
            return True
        elif self.f_hdf.isopen == 0:
            return False

    def open_file(self):
        """
        Checks if underlying snapshot file is open and available.  If not,
        opens the file.
        
        Returns
        --------
        f_hdf : ``PyTables file``
            HDF file handle.

        """
        if self.is_open() is True:
            logger.user_information('%s already open' % self.snap.fn_base)
            return self.f_hdf
        else:
            if os.path.exists(self.snap.filename):
                if valid_file(self.snap.filename):
                    self.f_hdf = pyT.openFile(self.snap.filename, 'a')
                else:
                    self.logger.error('File %s not HDF5 file.' %
                                                            self.snap.fn_base)
            else:
                self.logger.user_message('Creating %s' % self.snap.fn_base)
                self.f_hdf = pyT.openFile(self.snap.filename, 'a')
            # TODO: Do not return - need to abstract f_hdf out of Snap
            return self.f_hdf

    def hdf_create_array(self, dataArrayName, dataArray, dtype='float',
                                    location='Node', update_meta=True):
        """
        Create a new HDF array within the current dataset.

        Wrapping of ``PyTables.createArray(...)``

        Parameters
        ----------
        dataArrayName : string
            Name of newly created HDF array.
        dataArray : :class:`numpy.ndarray`
            Numpy array to be written to HDF file.
        dtype : string,  {'Int' | 'Float'}
            Data type.
        location : string,  {'Node' | 'Cell'}
            Location of scalar attribute within grid cell.
        update_meta : boolean,  [default: True]
            Flag for updating HDF metadata.
            (e.g. ``col_mapping``, ``vectors``, etc.)

        """
        dtype = dtype.lower()
        if 'float' in dtype:
            dtype = 'Float'
        elif 'int' in dtype:
            dtype = 'Int'
        else:
            self.logger.error("Invalid data type.  Must be 'int*' or 'float*'")

        if location in settings.LOCATIONS:
            location = location.title()
        else:
            # TODO: Must be 'Node' for atmowise dataset.
            self.logger.error("Data location must be either 'Node' or 'Cell'")

        # Enforcing dtype from function parameter
        dataArray = vectors.vector_type_conv(dataArray, str(dtype))

        self.snap.f_hdf.createArray(self.snap.active_dataset,
                                     str(dataArrayName), dataArray)
        if update_meta is True:
            self.hdf_update_col_mapping(dataArrayName, location=location,
                                                                dtype=dtype)

    def hdf_create_EArray(self, EArrayName, shapeTuple,
                               dtype='float', location='Node',
                               expectedrows=100000, update_meta=True):
        """
        Create new **extensible** HDF array within current dataset.

        Wrapping of ``PyTables.createEArray(...)``

        Parameters
        -----------
        EArrayName : string
            Name of extendable array to be created
        shapeTuple : iterable
            Shape of ``EArray``. Must have at least one value set to 0.
        dtype : string,  {'Int' | 'Float'}
            Data type string for array.
        location : string,  {'Node' | 'Cell'}
            Location of scalar attribute within grid cell.
        expectedrows : int,  [default=100000]
            Good guess at expected number of rows to improve speed.
        update_meta : boolean,  [default: True]
            Flag for updating HDF metadata.
            (e.g. ``col_mapping``, ``vectors``, etc.).

        """
        expectedrows = int(expectedrows)

        dtype = dtype.lower()
        if 'float' in dtype:
            atom = pyT.FloatAtom()
        elif 'int' in dtype:
            atom = pyT.IntAtom()
        else:
            self.logger.error("Invalid data type.  Must be 'int*' or 'float*'")

        if location in settings.LOCATIONS:
            location = location.title()
        else:
            # TODO: Must be 'Node' for atmowise dataset.
            self.logger.error("Data location must be either 'Node' or 'Cell'")

        self.snap.f_hdf.createEArray(self.snap.active_dataset,
                              str(EArrayName), atom, shapeTuple, expectedrows)

        if update_meta is True:
            self.hdf_update_col_mapping(EArrayName, dtype=dtype,
                                                  location=location)

    def hdf_EArray_append_data(self, EArrayName, dataArray):
        """
        Append new data from dataArray to existing extensible EArray.

        Parameters
        ----------
        EArrayName : string
            Name of extensible HDF ``EArray`` target.
        dataArray : :class:`numpy.ndarray`
            Data to be added.

        """
        #TODO: Add try/except to handle if target EArray doesn't yet exist.
        #TODO:  Catch difference in shape btw dataArray and EArray
        EArray = self.snap.f_hdf.getNode('/'+self.snap.active_dataset_name,
                                                                   EArrayName)

        dtype = str(EArray.atom).lower()
        if 'float' in dtype:
            dtype = 'float'
        if 'int' in dtype:
            dtype = 'int'

        # Enforcing dtype from function parameter
        dataArray = vectors.vector_type_conv(dataArray, str(dtype))
        EArray.append(dataArray)

    def hdf_create_dataset(self, new_name):
        """
        Create new dataset at the root of HDF file.

        Dataset, ``'rawSimResults'`` always exists to hold atomwise simulation
        results, so should not be used.

        .. note::
           Does not overwrite an existing dataset.

        Parameters
        ----------
        new_name : string
            Name of new data set to create at root (i.e. ``'/'``) of HDF file.

        """
        #TODO: Add try/except around attempt to create new dataset.
        self.snap.f_hdf.createGroup('/', str(new_name))

    def hdf_delete_array(self, arrayName):
        """
        Remove existing HDF array from current dataset.

        Parameters
        ----------
        arrayName : string
            Name of array to remove.

        """
        #TODO: Add try/except around attempt to get target array.
        array = self.snap.active_dataset._f_getChild(str(arrayName))
        array.remove()

    def hdf_rename_dataset(self, newName):
        """
        Renames currently active dataset on the fly so that snap
        continues to point to same data before and after.

        Parameters
        ----------
        newName : string
            New dataset name.

        """
        self.snap.active_dataset._f_rename(str(newName))
        self.flush()
        self.snap.set_active_dataset(str(newName))

    def hdf_delete_dataset(self, dataset_name):
        """
        Removes dataset located at the root of the HDF file.

        .. warning::
           Acts recursively and destructively, deleting all data arrays
           contained within the dataset.

        Parameters
        ----------
        dataset_name : string
            Name of dataset to delete.

        """
        if dataset_name == 'rawSimResults':
            self.logger.error('Cannot delete rawSimResults.')
        else:
            if self.verbose is True:
                if util.user_approve('Remove %s ?' % dataset_name) is True:
                    self.snap.f_hdf.removeNode('/'+str(dataset_name),
                                                                recursive=True)
                else:
                    self.logger.warning('%s not removed.' % dataset_name)

    def hdf_add_meta(self, metaName, metaValue):
        """
        Add metadata to HDF file.

        Parameters
        ----------
        metaName : string
            Name of metadata item.
        metaValue : object
            Object to store in HDF metadata for active dataset in HDF file.

        """
        setattr(self.snap.active_dataset._v_attrs, str(metaName), metaValue)

    def hdf_rename_meta(self, currentName, newName):
        """
        Rename existing metadata entry in HDF file.

        Parameters
        ----------
        currentName : string
            Current name of metadata entry.
        newName : string
            New name of metadata entry.

        """
        self.snap.active_dataset._v_attrs._f_rename(
                                                str(currentName), str(newName))

    def hdf_modify_meta(self, metaName, metaValue):
        """
        Change value of existing metadata entry in HDF file.

        Parameters
        ----------
        metaName : string
            Name of existing metadata entry.
        metaValue : object
            New value of metadata entry.

        """
        setattr(self.snap.active_dataset._v_attrs, str(metaName), metaValue)

    def hdf_delete_meta(self, metaName):
        """
        Delete existing metadata entry.

        Parameters
        ----------
        metaName : string
            Name of metadata entry to be deleted.

        """
        #TODO: try/except in case entry does not exist.
        self.snap.active_dataset._f_delAttr(metaName)

    def hdf_update_col_mapping(self, name, dtype='', location=''):
        """
        Update ``col_mapping`` metadata in current HDF dataset to reflect
        changes in underlying data arrays.

        Parameters
        ----------
        name : string
            Name of scalar attribute to update.
        dtype:  string,  {'Int' | 'Float'}, [default='']
            New dtype.  If not changing dtype, do not use argument.
        location : string,  {'Node' | 'Cell'},  [default='']
            New location of scalar attribute.
            If not changing location, do not use argument.

        """
        # TODO: checks for dtype and location w/i allowable options.
        if dtype == '':
            dtype = self.snap.col_mapping[name]['dtype']
        if location == '':
            location = self.snap.col_mapping[name]['location']

        curr_mapping = self.snap.col_mapping

        # Adds new entry to col_mapping if name is not present
        curr_mapping[name] = {'dtype': dtype, 'location': location}
        self.hdf_modify_meta('col_mapping', curr_mapping)

    def hdf_update_composite_mapping(self, name, type_mapping,
                                   dtype='', location='',  comps=''):
        """
        Updates the metadata used to describe the composition of
        vectors and tensors (e.g. components, dtype, grid location).

        Parameters
        ----------
        name : string
            Name of attribute to update (e.g. vector or tensor name).
        type_mapping : string,  {'vectors' |  'tensors' | 'symm_tensors'}
            Type of attribute.
        dtype : string,  {'Int' | 'Float'},  [default='']
            New dtype.  If not changing dtype, use default.
        location : string,  {'Node' | 'Cell'},  [default='']
            New location of scalar attribute within grid cell.
            If not changing location, use default.
        comps : iterable,  [default='']
            Name of underlying data used for components of vector or tensor.
            (e.g. ``['vx', 'vy', 'vz']``)
            Use ``None`` for a contiguous item without separate components
            (i.e. a vector stored as a (Nx3) array rather than 3 (Nx1) arrays).
            Use the default if want to leave ``comps`` unchanged.

        """
        if (type_mapping not in settings.PROP_SIZE):
            self.logger.error('Invalid collection type %s' % type_mapping)

        if type_mapping == 'scalars':
            self.logger.warning("Collection type is scalars - do nothing")
            return

        curr_mapping = self.snap.get_metadata(type_mapping)

        if dtype == '':
            dtype = curr_mapping[name]['dtype']
        if location == '':
            location = curr_mapping[name]['location']
        if comps == '':
            # Not a contiguous item i.e. separate comp names exist.
            if 'comps' in curr_mapping[name]:
                comps = curr_mapping[name]['comps']
            # Contiguous parameter w/o need for separate comp names.
            else:
                comps = None

        assert dtype in settings.DTYPES
        assert location in settings.LOCATION
        if comps is not None:
                assert len(comps) == settings.PROP_SIZE[type_mapping]

        curr_mapping[name] = {'location': location,
                              'dtype': dtype,
                              'comps': comps}

        self.hdf_modify_meta(type_mapping, curr_mapping)

    def flush(self):
        """
        Flusses I/O buffer to bring HDF file on disk up to date
        with any changes.

        Wrapper on PyTables flush function.

        """
        self.f_hdf.flush()

    def cleanup(self):
        """
        Cleanup and close out HDF file.

        .. note:: 
           When Python process exits, existing PyTables flush and close
           functions run automatically.

        """
        self.f_hdf.flush()
        self.f_hdf.close()
        self.logger.user_message('Flush and close %s' % self.snap.filename)




# =============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    main()
