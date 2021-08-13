from input2d import *
# Computing stiffness matrix with considering of distributed loads    
def AssemblePlaneBeam(K, Ke, LeftNode, RightNode):

    LeftDOFu = LeftNode * ndim
    LeftDOFv = LeftDOFu + 1
    LeftDOFtheta = LeftDOFu + 2
    
    RightDOFu = RightNode * ndim
    RightDOFv = RightDOFu + 1
    RightDOFtheta = RightDOFu + 2
    
    StiffnessGlobal = np.array([LeftDOFu, LeftDOFv, LeftDOFtheta, \
                                RightDOFu, RightDOFv, RightDOFtheta])
    print('Global element stiffness index array:', StiffnessGlobal)
    # Append new values into initial global stiffness matrix
    for i in StiffnessGlobal:
        for j in StiffnessGlobal:
            m = np.array(np.where(StiffnessGlobal==i)).flatten()[0]
            n = np.array(np.where(StiffnessGlobal==j)).flatten()[0]
            ##K[i,j] = K[i,j] + Ke[i-LeftDOFu,j-LeftDOFu]
            K[i,j] = K[i,j] + Ke[m,n]
    
    return K