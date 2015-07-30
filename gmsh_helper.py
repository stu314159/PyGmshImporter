"""
gmsh_helper.py

Written by: Stu Blair
Date: 6/17/2015

Purpose: Provide a module with various utilities needed for importing
GMSH mesh data into the python environment.

Last updated:


"""
from scipy.spatial import ConvexHull
import numpy.linalg as la
import numpy as np


class genElement:
    def __init__(self,num,elType,tags,nds):
        self.num = int(num)
        self.elType = int(elType)
        self.tags = tags
        self.nds = nds
        

class Node:
    def __init__(self,num,x,y,z):
        self.num = int(num)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

class Tet:
    def __init__(self,num,nd1,nd2,nd3,nd4):
        self.num = num
        self.nd1 = nd1
        self.nd2 = nd2
        self.nd3 = nd3
        self.nd4 = nd4
        self.D = np.zeros((4,4)).astype(np.float64)
        self.D[0][:] = [self.nd1.x, self.nd1.y, self.nd1.z, 1.]
        self.D[1][:] = [self.nd2.x, self.nd2.y, self.nd2.z, 1.]
        self.D[2][:] = [self.nd3.x, self.nd3.y, self.nd3.z, 1.]
        self.D[3][:] = [self.nd4.x, self.nd4.y, self.nd4.z, 1.]
        #self.sign = np.sign(la.det(self.D))
        self.D0 = la.det(self.D)
        self.sign = np.sign(self.D0)
        self.Ds = np.zeros((4,4,4)).astype(np.float64)
        for i in range(4):
            self.Ds[i][:][:] = np.copy(self.D)
        
        
    def pointInside(self,pt):
        """
        pt is a 3-tuple with the coordinates of a point that I would like to check to see if it is
        inside the tetrahedron or not.
        
        """
        
        for i in range(4):
            self.Ds[i][i][:3]=pt
            
        dets = np.zeros(4).astype(np.float64)
        sumD = 0.
        for i in range(4):
            dets[i]=la.det(self.Ds[i][:][:])
            if abs(dets[i]) < 1e-6: # prevent floating point roundoff from giving me a problem
                dets[i]=0.
            sumD+=dets[i]
            if dets[i]==0.:
                return True # point is coincident with a node and/or on a boundary
            if np.sign(dets[i]) == np.sign(-1*self.sign):
                return False #point is outside of 
        
        
        
        return True
        
        """
        This is the basic implementation of the convex hull idea.  It works, but is slow and is not amenable to
        acceleration via ocl tricks.
        
        #p1 = [self.nd1.x,self.nd1.y,self.nd1.z]
        #p2 = [self.nd2.x,self.nd2.y,self.nd2.z]
        #p3 = [self.nd3.x,self.nd3.y,self.nd3.z]
        #p4 = [self.nd4.x,self.nd4.y,self.nd4.z]
        #
        #conv_hull = ConvexHull([p1,p2,p3,p4],incremental=True)
        #v1 = conv_hull.vertices
        #
        #conv_hull.add_points([pt])
        #v2 = conv_hull.vertices
        #
        #return set(v1)==set(v2) # if the point is inside the tetrahedron, there will be no change in the set of points for the convex hull.
        """
        
def getTets(genElList,nodeList):
    """
    input: a list of elements (of type genElement)
    output: a list of tetrahedron elements (of type Tet)
    """
    tetList = []
    ind = 0
    for key in genElList:
        el = genElList[key]
        if el.elType == 4:
            nd1 = nodeList[el.nds[0]]
            nd2 = nodeList[el.nds[1]]
            nd3 = nodeList[el.nds[2]]
            nd4 = nodeList[el.nds[3]]
            tetList.append(Tet(ind,nd1,nd2,nd3,nd4)) #<-- this is a truly crappy design.
            ind+=1
    return tetList
        
        

def msh_import(msh_filename):
    """
    msh_filename = string indicating the name of the *.msh file
    
    
    """
   
    with open(msh_filename) as msh_file:
        # skip text before begining of interesting block
        for line in msh_file:
            if line.strip() == '$Nodes':
                break
        #numNd = msh_file.readline()
        #assert(len(numNd)==1)
        #numNd = int(numNd)
        #print 'There are %d nodes to read.'%numNd
        nodes = {}
        for line in msh_file:
            if line.strip()=='$EndNodes':
                break
            node_entries = line.split()
            if len(node_entries)==1:
                numNd=int(node_entries[0])
                print 'There are %d nodes to read.'%numNd
                continue
            
            nodeNum = int(node_entries[0])
            nodes[nodeNum]=Node(nodeNum,node_entries[1],node_entries[2],node_entries[3])
            #nodeList.append(Node(int(node_entries[0]),float(node_entries[1]),\
            #                         float(node_entries[2]),float(node_entries[3])))
        
        
        
        for line in msh_file:
            if line.strip() == '$Elements':
                break
        
                
        genEl = {}
        for line in msh_file:
            if line.strip() == '$EndElements':
                break
            
            elDat = line.split()
            
            if len(elDat)==1:
                numEl = int(elDat[0])
                print 'There are %d elements to read.'%numEl
                continue
            
            elNum = int(elDat[0])
            elType = int(elDat[1])
            elNumTags = int(elDat[2])
            elTags = []
            for i in range(3,3+elNumTags):
                elTags.append(int(elDat[i]))
                
            # remaining elDat are node numbers
            nds = []
            ndList = elDat[(3+elNumTags):]
            for nd in ndList:
                nds.append(int(nd))
            
            genEl[elNum]=genElement(elNum,elType,elTags,nds)
            #genEl.append(genElement(elNum,elType,elTags,nds))
    
    return nodes,genEl
    
    
    
    
if __name__=="__main__":
  nodeList,genElementList= msh_import("WallMountedBrickSW.msh")
  tetList = getTets(genElementList,nodeList)
  print 'All Done!!'