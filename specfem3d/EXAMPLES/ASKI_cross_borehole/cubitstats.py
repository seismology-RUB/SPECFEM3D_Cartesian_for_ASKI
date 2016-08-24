#----------------------------------------------------------------------------
#   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.0.
#
#   ASKI version 1.0 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.0 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.0.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#
import cubit
from math import sqrt
from numpy import histogram
############################# FUNCTIONS #########################
class cubitstats():
    def __init__(self):
        print('pass')
        self.getstats()

    def minEdgeHex(self,x):
        ''' Take Element and return min edge length.
        x has to be a matrix of eight xyz coordinates'''
        length=[]
        length.append(sqrt((x[0][0]-x[1][0])**2 + (x[0][1]-x[1][1])**2 + (x[0][2]-x[1][2])**2 ))
        length.append(sqrt((x[1][0]-x[2][0])**2 + (x[1][1]-x[2][1])**2 + (x[1][2]-x[2][2])**2 ))
        length.append(sqrt((x[2][0]-x[3][0])**2 + (x[2][1]-x[3][1])**2 + (x[2][2]-x[3][2])**2 ))
        length.append(sqrt((x[3][0]-x[0][0])**2 + (x[3][1]-x[0][1])**2 + (x[3][2]-x[0][2])**2 ))
        length.append(sqrt((x[4][0]-x[5][0])**2 + (x[4][1]-x[5][1])**2 + (x[4][2]-x[5][2])**2 ))
        length.append(sqrt((x[5][0]-x[6][0])**2 + (x[5][1]-x[6][1])**2 + (x[5][2]-x[6][2])**2 ))
        length.append(sqrt((x[6][0]-x[7][0])**2 + (x[6][1]-x[7][1])**2 + (x[6][2]-x[7][2])**2 ))
        length.append(sqrt((x[7][0]-x[4][0])**2 + (x[7][1]-x[4][1])**2 + (x[7][2]-x[4][2])**2 ))
        length.append(sqrt((x[0][0]-x[4][0])**2 + (x[0][1]-x[4][1])**2 + (x[0][2]-x[4][2])**2 ))
        length.append(sqrt((x[1][0]-x[5][0])**2 + (x[1][1]-x[5][1])**2 + (x[1][2]-x[5][2])**2 ))
        length.append(sqrt((x[2][0]-x[6][0])**2 + (x[2][1]-x[6][1])**2 + (x[2][2]-x[6][2])**2 ))
        length.append(sqrt((x[3][0]-x[7][0])**2 + (x[3][1]-x[7][1])**2 + (x[3][2]-x[7][2])**2 ))
        l=sorted(length)
        return l[0]

    def maxEdgeHex(self,x):
        ''' Take Element and return may edge length.
        x has to be a matrix of eight xyz coordinates'''
        length=[]
        length.append(sqrt((x[0][0]-x[1][0])**2 + (x[0][1]-x[1][1])**2 + (x[0][2]-x[1][2])**2 ))
        length.append(sqrt((x[1][0]-x[2][0])**2 + (x[1][1]-x[2][1])**2 + (x[1][2]-x[2][2])**2 ))
        length.append(sqrt((x[2][0]-x[3][0])**2 + (x[2][1]-x[3][1])**2 + (x[2][2]-x[3][2])**2 ))
        length.append(sqrt((x[3][0]-x[0][0])**2 + (x[3][1]-x[0][1])**2 + (x[3][2]-x[0][2])**2 ))
        length.append(sqrt((x[4][0]-x[5][0])**2 + (x[4][1]-x[5][1])**2 + (x[4][2]-x[5][2])**2 ))
        length.append(sqrt((x[5][0]-x[6][0])**2 + (x[5][1]-x[6][1])**2 + (x[5][2]-x[6][2])**2 ))
        length.append(sqrt((x[6][0]-x[7][0])**2 + (x[6][1]-x[7][1])**2 + (x[6][2]-x[7][2])**2 ))
        length.append(sqrt((x[7][0]-x[4][0])**2 + (x[7][1]-x[4][1])**2 + (x[7][2]-x[4][2])**2 ))
        length.append(sqrt((x[0][0]-x[4][0])**2 + (x[0][1]-x[4][1])**2 + (x[0][2]-x[4][2])**2 ))
        length.append(sqrt((x[1][0]-x[5][0])**2 + (x[1][1]-x[5][1])**2 + (x[1][2]-x[5][2])**2 ))
        length.append(sqrt((x[2][0]-x[6][0])**2 + (x[2][1]-x[6][1])**2 + (x[2][2]-x[6][2])**2 ))
        length.append(sqrt((x[3][0]-x[7][0])**2 + (x[3][1]-x[7][1])**2 + (x[3][2]-x[7][2])**2 ))
        l=sorted(length)
        return l[11]

    def calcMaxFreq(self,x,vs):
        ''' For a given edge length x and a given velocity vs, get the max frequency resolved by the mesh'''
        return vs/(x*1.25) # compare SPECFEM3D-3.0/src/shared/check_mesh_resolution.f90 lines 215,204: 1.25 = 5/4 , where 5 is number of GLL points per wavelength and 4 = NGLL-1 is the number of intervals between GLL points in an element

    def calcMaxDt(self,x,vp,cn):
        ''' Calc the max dt for the given Block, with x= max_edgelengh, vp'''
        dt=x*0.173*cn/vp
        return dt

    def getstats(self):
        ''' Get statistics for differend blocks in the hexmesh '''
    ############################## READ BLOCKS ######################
        try:
            blocks=cubit.get_block_id_list()
        except:
            print('no blocks defined')

        nblocks=len(blocks)
        print('Found ',nblocks,' Blocks')
        # lists of details of all elastic blocks. For each elastic block, append respective values
        name=[]
        typ=[]
        nattrib=[]
        vp=[]
        vs=[]
        rho=[]
        minEdge=[]
        maxEdge=[]
        minEdgeMesh = []
        maxEdgeMesh = []
        # k: counter of all elastic blocks
        k=0
        for block in blocks:
            minEdgeMeshTmp = []
            maxEdgeMeshTmp = []
            blockname=cubit.get_exodus_entity_name('block',block)
            # only apply stats to elastic blocks
            if blockname.find("elastic") >= 0:
                # NOW APPEND VALUES OF THIS ELASTIC BLOCK TO LISTS (defined as empty above)
                name.append(blockname)
                typ.append(cubit.get_block_element_type(block))
                # attribute count
                nattrib.append(cubit.get_block_attribute_count(block))
                # vp
                vp.append(cubit.get_block_attribute_value(block,1))
                # vs
                vs.append(cubit.get_block_attribute_value(block,2))
                # rho
                rho.append(cubit.get_block_attribute_value(block,3))
                hexes=cubit.get_block_hexes(block)
                print(block, k,'name: ',name[k],'#hexes:',len(hexes),'type: ',typ[k],'attribute count: ',nattrib[k],'vp: ',vp[k],'vs: ',vs[k],'rho: ',rho[k])
                
                # minimum/maximum edgelength search for this elastic block
                # initiate  with very large minEdge and very small maxEdge
                minEdge.append(1e9)
                maxEdge.append(1e-9)
                # now search
                for hexa in hexes:
                    nodes=cubit.get_connectivity('Hex',hexa)
                    # get coordinates of nodes of this hexahedron
                    node_coords=[]
                    for node in nodes:
                        node_coords.append(cubit.get_nodal_coordinates(node))

                    minv=self.minEdgeHex(node_coords)
                    minEdgeMeshTmp.append(minv)
                    maxv=self.maxEdgeHex(node_coords)
                    maxEdgeMeshTmp.append(maxv)

                    # is minimum edgelength of this this element smaller than minimum edgelength found before?
                    if minv < minEdge[k]:
                        minEdge[k]=minv

                    # is maximum edgelength of this this element larger than maximum edgelength found before?
                    if maxv > maxEdge[k]:
                        maxEdge[k]=maxv

                minEdgeMesh.append(minEdgeMeshTmp)
                maxEdgeMesh.append(maxEdgeMeshTmp)
                # end of processing of this elastic block. incremet counter
                k=k+1
            # if not "elastic" in this blockname
            else:
                print(blockname, '--no elastic')

        # now calculate maximum frequency for which the mesh is stable
        # and maximal dt

        # first do that for each block
        fmax=[]
        dtmax_05=[]
        for i in range(len(name)):
            fmax.append(self.calcMaxFreq(maxEdge[i],vs[i]))
            dtmax_05.append(self.calcMaxDt(minEdge[i],vp[i],cn=0.5))
            min_dx = minEdge[i]*0.1727
            max_dx = maxEdge[i]*0.3272*sqrt(3.)
            print('Block ', name[i], ' min dx: ',min_dx, ' max dx: ',max_dx, 'max frequency: ',fmax[i],
                  ' maximal dt (C=0.5):',dtmax_05[i])
            vals,bin_edges = histogram(minEdgeMesh[i],10)
            nhexa = sum(vals)
            print('   '+str(nhexa)+' HEXAS IN THIS BLOCK')
            print('   histogram of minimum edge length of hexa of this block:')
            for j,val in enumerate(vals):
                print('       %6.3f %% of hexas have minimum edgelength between %3.3f and %3.3f  ->   dt (C=0.5) is %8.3e'%(100.*val/float(nhexa),bin_edges[j],bin_edges[j+1],self.calcMaxDt(bin_edges[j],vp[i],cn=0.5)))
            vals,bin_edges = histogram(maxEdgeMesh[i],10)
            nhexa = sum(vals)
            print('   '+str(nhexa)+' HEXAS IN THIS BLOCK')
            print('   histogram of maximum edge length of hexa of this block:')
            for j,val in enumerate(vals):
                print('       %6.3f %% of hexas have maximum edgelength between %3.3f and %3.3f  ->   fmax (5 GLL points per wavelength) is %8.3e'%(100.*val/float(nhexa),bin_edges[j],bin_edges[j+1],self.calcMaxFreq(bin_edges[j+1],vs[i])))
        # now compare the values for all the blocks
        dtm_05=sorted(dtmax_05)
        print('The minimum over all blocks of the respective maximal dt for C=0.5 is: ',dtm_05[0])
        fmax_min=sorted(fmax)
        print('The minimum over all blocks of the respective maximal frequency is: ',fmax_min[0])

def stats():
    x=cubitstats()

def getEdgeLength(fmax,vmin):
    return vmin/(fmax*1.25) # compare SPECFEM3D-3.0/src/shared/check_mesh_resolution.f90 lines 215,204: 1.25 = 5/4 , where 5 is number of GLL points per wavelength and 4 = NGLL-1 is the number of intervals between GLL points in an element

def getNSAMP(T,dt):
    return T/dt
        
if __name__ == '__main__':
    stats()
    




