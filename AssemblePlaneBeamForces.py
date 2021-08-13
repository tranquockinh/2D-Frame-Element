from input2d import *
#Assemble the prescribed concentrated loads
def AssemblePlaneBeamForces(p, cloads):

    # p is a vector size (ndof,1)
    # cloads is a matrix defined with loading components on loaded 
    #node numbers
    
    numLoadedNodes = np.size(cloads, 0) ## number of nodes loaded
    
    for i in range(numLoadedNodes):
        ## Pull load components of each current node
        nodenum = cloads[i, 0]
        px = cloads[i, 1]
        py = cloads[i, 2]
        mz = cloads[i, 3]
        pNode = np.array([px, py, mz])
        ## pass to global concentrated load vector
        ## Locate the coordinates of global loading components
        pxDof = nodenum * ndim
        pyDof = nodenum * ndim + 1
        mzDof = nodenum * ndim + 2
        pDof = np.array([pxDof, pyDof, mzDof], dtype=int)
        ##print('pDof is:\n', pDof)
        ## Append the load components of each node to the right coordinates
        for j in pDof:
            p[j] = p[j] + pNode[j-int(pxDof)]
    
    return p
