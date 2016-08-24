#!/usr/bin/env python
#!python

import cubit
import cubit2specfem3d
import cubitstats

reload(cubit2specfem3d)
reload(cubitstats)

import os
import sys

cubit.cmd('#####################################################')
cubit.init([""])
cubit.cmd('reset')


############################
cubit.cmd('## CREATE BASE BRICK')

cubit.cmd('brick x 200 y 200 z 150')
cubit.cmd('volume 1 move x 0.0 y 0.0 z -75.0')

#cubit.cmd('imprint all')
#cubit.cmd('merge all')

############################
cubit.cmd('## START MESHING')

# REQUESTED STABILITY:
# indicate final minimal velocity in the model
# indicate maximum stable frequency
vmin = 625.0
fmax = 36.0
# COMPUTE ELEMENT SIZE:
#el = cubitstats.getEdgeLength(fmax,vmin)
el = 4.0 # for tenth ASKI iteration, same as in dissertation, and gji aski paper. this assures a good coverage of wavefield points

print ("## el = ",el)
cubit.cmd('vol all size %10.10f'%el)

cubit.cmd('mesh vol all')

cubit.cmd('## END OF MESHING')
############################


############################
cubit.cmd('## ASSIGNING MATERIAL PROPERTIES')
#
# HALFSPACE
cubit.cmd('block 1 hex in vol 1')
cubit.cmd('block 1 name "elastic 1"         # Material 1, background')
cubit.cmd('block 1 attribute count 7')
cubit.cmd('block 1 attribute index 1 1     # flag for material: 1 for 1. material')
cubit.cmd('block 1 attribute index 2 1082.0  # vp')
cubit.cmd('block 1 attribute index 3 625.0  # vs')
cubit.cmd('block 1 attribute index 4 1800.0  # rho')
cubit.cmd('block 1 attribute index 5 9999.0')  # Qkappa
cubit.cmd('block 1 attribute index 6 9999.0')  # Qmu
cubit.cmd('block 1 attribute index 7 0       # anisotropy_flag')
#
#
cubit.cmd('## BOUNDARY CONDITIONS')
cubit.cmd('## free surface')
cubit.cmd('block 2 face in surface 1')
cubit.cmd('block 2 name "face_topo"')

cubit.cmd('## absorbing bottom')
cubit.cmd('block 3 face in surface 2')
cubit.cmd('block 3 name "face_abs_bottom"')

cubit.cmd('## absorbing ymin:')
cubit.cmd('block 4 face in surface 3')
cubit.cmd('block 4 name "face_abs_ymin"')

cubit.cmd('## absorbing xmin:')
cubit.cmd('block 5 face in surface 4')
cubit.cmd('block 5 name "face_abs_xmin"')

cubit.cmd('## absorbing ymax:')
cubit.cmd('block 6 face in surface 5')
cubit.cmd('block 6 name "face_abs_ymax"')

cubit.cmd('## absorbing xmax:')
cubit.cmd('block 7 face in surface 6')
cubit.cmd('block 7 name "face_abs_xmax"')

cubitstats.stats()

os.system('mkdir -p MESH')

cubit.cmd('save as "MESH/meshing.cub" overwrite')

cubit2specfem3d.export2SESAME('MESH')





