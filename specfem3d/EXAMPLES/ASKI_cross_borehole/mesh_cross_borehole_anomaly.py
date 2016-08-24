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


############################
cubit.cmd('## CREATE ANOMALY BRICK')

#cubit.cmd('## WEBCUTS X-DIRECTION FOR ANOMALY')
cubit.cmd('webcut vol all with plane xplane offset -15.0 noimprint nomerge')
cubit.cmd('webcut vol all with plane xplane offset  5.0 noimprint nomerge')
cubit.cmd('webcut vol all with plane xplane offset 15.0 noimprint nomerge')

#cubit.cmd('## WEBCUTS Y-DIRECTION FOR ANOMALY')
cubit.cmd('webcut vol all with plane yplane offset -15.0 noimprint nomerge')
cubit.cmd('webcut vol all with plane yplane offset 5.0 noimprint nomerge')
cubit.cmd('webcut vol all with plane yplane offset 15.0 noimprint nomerge')

#cubit.cmd('## WEBCUTS Z-DIRECTION FOR ANOMALY')
cubit.cmd('webcut vol all with plane zplane offset -40.0 noimprint nomerge')
cubit.cmd('webcut vol all with plane zplane offset -25.0 noimprint nomerge')
cubit.cmd('webcut vol all with plane zplane offset -10.0 noimprint nomerge')

cubit.cmd('imprint all')
cubit.cmd('merge all')

############################
cubit.cmd('## START MESHING')

# REQUESTED STABILITY:
# indicate final minimal velocity in the model
# indicate maximum stable frequency
vmin = 575.0
fmax = 150.0 # ricker with f_c = 50 Hz , request stability of 3*f_c = 150 Hz
# COMPUTE ELEMENT SIZE:
el = cubitstats.getEdgeLength(fmax,vmin)

print ("## el = ",el)
cubit.cmd('vol all size %10.10f'%el)

cubit.cmd('mesh vol all')

cubit.cmd('## END OF MESHING')
############################


############################
cubit.cmd('## ASSIGNING MATERIAL PROPERTIES')
#
# BACKGROUND IS EVERYTHING EXCEPT VOL 39 40 43 44 55 56 59 60
cubit.cmd('block 1 hex in vol 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 45 46 47 48 49 50 51 52 53 54 57 58 61 62 63 64')
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
# ANOMALY BOX 1, - fast block - 108% of background
cubit.cmd('block 2 hex in vol 43 44 59 60')
cubit.cmd('block 2 name "elastic 2"         # Material 2')
cubit.cmd('block 2 attribute count 7')
cubit.cmd('block 2 attribute index 1 2     # flag for material: 1 for 1. material')
cubit.cmd('block 2 attribute index 2 1169.0  # vp')
cubit.cmd('block 2 attribute index 3 675.0  # vs')
cubit.cmd('block 2 attribute index 4 1800.0  # rho')
cubit.cmd('block 2 attribute index 5 9999.0')  # Qkappa
cubit.cmd('block 2 attribute index 6 9999.0')  # Qmu
cubit.cmd('block 2 attribute index 7 0       # anisotropy_flag')
#
# ANOMALY BOX 2, - medium fast block - 105% of background
cubit.cmd('block 3 hex in vol 39 40')
cubit.cmd('block 3 name "elastic 3"         # Material 3')
cubit.cmd('block 3 attribute count 7')
cubit.cmd('block 3 attribute index 1 3     # flag for material: 1 for 1. material')
cubit.cmd('block 3 attribute index 2 1136.0  # vp')
cubit.cmd('block 3 attribute index 3 656.0  # vs')
cubit.cmd('block 3 attribute index 4 1800.0  # rho')
cubit.cmd('block 3 attribute index 5 9999.0')  # Qkappa
cubit.cmd('block 3 attribute index 6 9999.0')  # Qmu
cubit.cmd('block 3 attribute index 7 0       # anisotropy_flag')
#
# ANOMALY BOX 3, - medium slow block - 95% of background
cubit.cmd('block 4 hex in vol 55')
cubit.cmd('block 4 name "elastic 4"         # Material 4')
cubit.cmd('block 4 attribute count 7')
cubit.cmd('block 4 attribute index 1 4     # flag for material: 1 for 1. material')
cubit.cmd('block 4 attribute index 2 1028.0  # vp')
cubit.cmd('block 4 attribute index 3 594.0  # vs')
cubit.cmd('block 4 attribute index 4 1800.0  # rho')
cubit.cmd('block 4 attribute index 5 9999.0')  # Qkappa
cubit.cmd('block 4 attribute index 6 9999.0')  # Qmu
cubit.cmd('block 4 attribute index 7 0       # anisotropy_flag')
#
# ANOMALY BOX 4, - slow block - 92% of background
cubit.cmd('block 5 hex in vol 56')
cubit.cmd('block 5 name "elastic 5"         # Material 5')
cubit.cmd('block 5 attribute count 7')
cubit.cmd('block 5 attribute index 1 5     # flag for material: 1 for 1. material')
cubit.cmd('block 5 attribute index 2 995.0  # vp')
cubit.cmd('block 5 attribute index 3 575.0  # vs')
cubit.cmd('block 5 attribute index 4 1800.0  # rho')
cubit.cmd('block 5 attribute index 5 9999.0')  # Qkappa
cubit.cmd('block 5 attribute index 6 9999.0')  # Qmu
cubit.cmd('block 5 attribute index 7 0       # anisotropy_flag')


cubit.cmd('## BOUNDARY CONDITIONS')
cubit.cmd('## free surface')
cubit.cmd('block 6 face in surface 51 61 71 39 90 100 110 78 129 139 149 121 135 145 155 123 ')
cubit.cmd('block 6 name "face_topo"')

cubit.cmd('## absorbing bottom')
cubit.cmd('block 7 face in surface 133 143 153 125 131 141 151 119 88 98 108 80 49 59 69 41 ')
cubit.cmd('block 7 name "face_abs_bottom"')

cubit.cmd('## absorbing ymin:')
cubit.cmd('block 8 face in surface 176 186 196 166 334 344 354 324 495 505 515 485 488 498 508 478 ')
cubit.cmd('block 8 name "face_abs_ymin"')

cubit.cmd('## absorbing xmin:')
cubit.cmd('block 9 face in surface 295 253 213 173 453 416 376 336 616 573 533 493 609 571 531 491 ')
cubit.cmd('block 9 name "face_abs_xmin"')

cubit.cmd('## absorbing ymax:')
cubit.cmd('block 10 face in surface 286 316 306 296 444 474 464 454 605 635 625 615 598 628 618 608 ')
cubit.cmd('block 10 name "face_abs_ymax"')

cubit.cmd('## absorbing xmax:')
cubit.cmd('block 11 face in surface 165 206 246 283 323 363 403 446 486 526 566 603 479 519 559 601 ')
cubit.cmd('block 11 name "face_abs_xmax"')

cubitstats.stats()

os.system('mkdir -p MESH')

cubit.cmd('save as "MESH/meshing.cub" overwrite')

cubit2specfem3d.export2SESAME('MESH')







