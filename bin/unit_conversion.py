#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Quickly loads the unit conversion objects within the `units` module.

.. warning::
   Default values are used for :class:`~mdhandle.units.LJUnits` and
   :class:`~mdhandle.units.MetalUnits`.

"""

from mdhandle import units

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    lj_2_si = units.LJUnits()
    meta_2_lj = units.MetalUnits()
