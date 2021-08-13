import numpy as np

# Defind the structure coordinates

ndcoords = 12 * np.matrix([[0.0, 0.0],
                           [0.0, 10],
                           [10, 10],
                           [10, 0.0]])
##print ('Coordinates: \n', ndcoords) 
# Defining nodal connection
elems = np.matrix([[0, 1],
                   [1, 2],
                   [2, 3]])
##print ('Element connection: \n', elems) # checking
# Material properties
eparray = 30E6 * np.matrix([[10, 200],
                            [10, 100],
                            [10, 200]])
##print ('Element properties: \n', eparray) # checking

# Initialization
ndim = 3
numnodes = np.size(ndcoords, 0) # number of elements in the first column
ndof = ndim * numnodes
nelems = np.size(elems, 0)
print ('number of element is \n ', nelems)
print ('Number of nodes is \n', numnodes) 

K = np.zeros((ndof, ndof),dtype=float)
ke = np.zeros((ndim * 2, ndim * 2),dtype=float)
##print ('K is \n {} \n and \n ke is \n {}'.format(K, ke)) # Checking

p = np.zeros((ndof, 1)) # to ensure the force vector equals ndof
peq = np.zeros((ndof, 1))
pe = np.zeros((ndim * 2, 1))
u = np.zeros((ndof,1))
##print (p,'\n',peq,'\n',pe,'\n','\n',u)

# Boundary conditions
bc = np.matrix([[0, 0, 0],
                [0, 1, 0],
                [0, 2, 0],
                [3, 0, 0],
                [3, 1, 0],
                [3, 2, 0]], dtype=int)
print (bc)                
cloads = np.matrix([[1,     10000,      0,      0],
                    [2,     0,          0,      5000]])

print (cloads)

# Distributed load

qloads = np.matrix([[1,      0,      0,      0],
                    [2,      0,      0,      0]])


# Computing stiffness matrix without considering of distributed loads
def StiffnessPlaneBeamBasic(ex, ey, ep):
    EA = ep[0, 0]
    EI = ep[0, 1]
    b = np.matrix([[ex[1] - ex[0]],
                   [ey[1] - ey[0]]])
    
    L = float(np.sqrt(b.T * b))
    n = (1 / L) * b
    
    
    # Compute the transfomation matrix T for coordinate transformation 
    #in 2D space 
    T = np.matrix([[ n[0,0],     n[1,0],    0,      0,         0,        0],
                   [-n[1,0],     n[0,0],    0,      0,         0,        0],
                   [ 0,          0,         1,      0,         0,        0],
                   [ 0,          0,         0,      n[0,0],    n[1,0],   0],
                   [ 0,          0,         0,     -n[1,0],    n[0,0],   0],
                   [ 0,          0,         0,      0,         0,        1]])
    # Local element stiffness matrix initialization
    kle = np.zeros((ndof, ndof),dtype=float)
    kaa = float(EA / L)
    kvv = float(12 * EI / L**3)
    kvm = float(6 * EI / L**2)
    kmm = float(4 * EI / L)
    kmv = float(2 * EI / L)
    # Local element stiffness matrix
    kle = np.matrix([[kaa,    0,    0,    -kaa,   0,    0  ],
                     [ 0,     kvv,  kvm,   0,    -kvv,  kvm],
                     [ 0,     kvm,  kmm,   0,    -kvm,  kmv],
                     [-kaa,   0,    0,     kaa,   0,    0  ],
                     [ 0,    -kvv, -kvm,   0,     kvv, -kvm],
                     [ 0,     kvm,  kmv,   0,    -kvm,  kmm]])
    
    Ke = np.transpose(T) * kle * T
    
    return Ke, T, kle, L, b
    
# Computing stiffness matrix with considering of distributed loads    
def StiffnessPlaneBeam(ex, ey, ep, qloads):
    EA = ep[0, 0]
    EI = ep[0, 1]
    b = np.matrix([[ex[1] - ex[0]],
                   [ey[1] - ey[0]]])
    
    L = float(np.sqrt(b.T * b))
    n = (1 / L) * b
    
    # Compute the transfomation matrix T for coordinate transformation 
    #in 2D space 
    T = np.matrix([[ n[0,0],     n[1,0],    0,      0,         0,        0],
                   [-n[1,0],     n[0,0],    0,      0,         0,        0],
                   [ 0,          0,         1,      0,         0,        0],
                   [ 0,          0,         0,      n[0,0],    n[1,0],   0],
                   [ 0,          0,         0,     -n[1,0],    n[0,0],   0],
                   [ 0,          0,         0,      0,         0,        1]])
    
    lload = np.matrix([[qloads[0,0]],[qloads[0,1]],[qloads[0,2]]])
    qxle = float(T[0, :3] * lload)
    qyle = float(T[1, :3] * lload)
    mzle = float(T[2, :3] * lload)
    ## Compute local distributed load matrix
    fle = np.matrix([[qxle * L / 2],
                     [qyle * L / 2],
                     [qyle * L**2 / 12 - mzle * L / 2],
                     [qxle * L / 2],
                     [qyle * L / 2],
                     [-qyle * L**2 / 12 + mzle * L / 2]])
                     
                     ## Compute the global load vector
    pe = np.transpose(T) * fle
                   
    # Local element stiffness matrix initialization
    kle = np.zeros((ndof, ndof),dtype=float)
    kaa = float(EA / L)
    kvv = float(12 * EI / L**3)
    kvm = float(6 * EI / L**2)
    kmm = float(4 * EI / L)
    kmv = float(2 * EI / L)
    # Local element stiffness matrix
    kle = np.matrix([[kaa,    0,    0,    -kaa,   0,    0  ],
                     [ 0,     kvv,  kvm,   0,    -kvv,  kvm],
                     [ 0,     kvm,  kmm,   0,    -kvm,  kmv],
                     [-kaa,   0,    0,     kaa,   0,    0  ],
                     [ 0,    -kvv, -kvm,   0,     kvv, -kvm],
                     [ 0,     kvm,  kmv,   0,    -kvm,  kmm]])
    
    Ke = np.transpose(T) * kle * T
    
    return Ke, T, kle, pe, L, b   

    
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
p = AssemblePlaneBeamForces(p, cloads)
print('Global force vector is p = \n {}'.format(p))    

# Assemble stiffness matrix
for ThisElement in range(nelems):
    LeftNode = elems[ThisElement, 0]
    RightNode = elems[ThisElement, 1]
    ex = [ndcoords[LeftNode, 0], ndcoords[RightNode, 0]]
    ey = [ndcoords[LeftNode, 1], ndcoords[RightNode, 1]]
    ep = eparray[ThisElement, :]
    ##print('==========================', ep)
    
    ## Take care of distributed loads
    found = 0 ## suppose there are no distributed loads applied
    for i in range(np.size(qloads,0)):
        if ThisElement == qloads[i,0]:
            found = 1
            break
    if found:
    
        Ke, T, kle, pe, L, b = StiffnessPlaneBeam(ex, ey, ep, qloads[i, 1:])
        
        pe = np.array(pe.T).flatten() ## Fattenning pe to an array
        peLeft = np.append([LeftNode], pe[:3])
        peRight = np.append([RightNode], pe[3:])
        eqvcloads = np.matrix([peLeft, peRight])
        print('Equivalent loads:\n', eqvcloads)
        print(eqvcloads.dtype)
        ## Equivalent loads
        peq = AssemblePlaneBeamForces(peq, eqvcloads)
    else:        
        ## Call the function of stiffness matrix
        Ke,T, kle, L, b = StiffnessPlaneBeamBasic(ex, ey, ep) 
        
        ## Check symmetry of the stiffness matrix
        for i in range(2 * ndim):
            for j in range(2 * ndim):
                roundoff_error = 1.0E-12
                if i != j:
                    if (Ke[i, j] - Ke[j, i]) > roundoff_error:                    
                        print ('Check symmetric matrix')
                    else: continue
    ## Assemble plane beam
    K = AssemblePlaneBeam(K, Ke, LeftNode, RightNode)
    for i in range(ndof):
            for j in range(ndof):
                roundoff_error = 1.0E-12
                if i != j:
                    if (K[i, j] - K[j, i]) > roundoff_error:                    
                        print ('Check symmetric matrix')
                    else: continue
    ## Check local striffness matrix
    print(' \
    Global element stiffness matrix of element Ke{} =\n {} \n \
    transformation matrix T{} =\n {}'.format(ThisElement,Ke,ThisElement,T)
    ) ## Check global stiffness matrix
    print('Global stiffness matrix of element {} is K =\n {}'. \
    format(ThisElement, K))
    print('Shape of K is: \n', np.shape(K))
    print('length of element L{} is\n {}\n and b{} is\n {} '. \
    format(ThisElement,L,ThisElement,b))
    print(ex,ey,ep)
    print()
    ##print('Ke{}=\n{}'.format(ThisElement, np.transpose(T) * kle * T))
 
# Solving the displacement and global force vector
def Solver(K, u, p, peq, bc):
    numPrescribedDof = np.size(bc, 0)
    numFreeDof = ndof - numPrescribedDof
    idxPrescribedDof = np.zeros((numPrescribedDof, 1), dtype=int)
    idxFreeDof = np.zeros((numFreeDof, 1), dtype=int)
    
    ## Index of prescribed DOF
    for i in range(numPrescribedDof):
        nodePrescibed = bc[i, 0]
        xORyORtheta = bc[i, 1]
        NumericalValue = bc[i, 2]
        idxPrescribed = nodePrescibed * ndim + xORyORtheta
        u[idxPrescribed] = NumericalValue
        idxPrescribedDof[i,0] = idxPrescribed
    
    
    ## Index of free DOF
    FreeDof = np.arange(0,ndof)
    FreeDof = np.delete(FreeDof, idxPrescribedDof)
    
    ## Compute global displacement and global nodal forces
    
    ## Initialization
    
    u_compute = np.zeros((np.size(FreeDof),1), dtype=float)
    K_compute = np.zeros((np.size(FreeDof),np.size(FreeDof)))
    K_prescribed = np.zeros((np.size(idxPrescribedDof),np.size(idxPrescribedDof)))
    K_Freepres = np.zeros((np.size(FreeDof),np.size(idxPrescribedDof)))
    K_presFree = np.zeros((np.size(idxPrescribedDof),np.size(FreeDof)))
    p_compute = np.zeros((np.size(FreeDof),1), dtype=float)
    peq_compute = np.zeros((np.size(FreeDof),1), dtype=float)
    peq_prescribed = np.zeros((np.size(idxPrescribedDof),1), dtype=float)
    u_Prescribed = np.zeros((np.size(idxPrescribedDof),1), dtype=float)           
    
    for j in FreeDof:
        for k in FreeDof:   
            m = np.array(np.where(FreeDof == j)).flatten()[0]
            n = np.array(np.where(FreeDof == k)).flatten()[0]
            p_compute[m] = p[j]
            peq_compute[m] = peq[j]
            K_compute[m, n] = K[j, k]
    
    for j in idxPrescribedDof:
        for k in idxPrescribedDof:
            m = np.array(np.where(idxPrescribedDof == j)).flatten()[0]
            n = np.array(np.where(idxPrescribedDof == k)).flatten()[0]
            peq_prescribed[m] = p[j]
            u_Prescribed[m] = u[j]
            K_prescribed[m,n] = K[j,k]
            
    for j in idxPrescribedDof:
        for k in FreeDof:
            m = np.array(np.where(idxPrescribedDof == j)).flatten()[0]
            n = np.array(np.where(FreeDof == k)).flatten()[0]
            K_presFree[m,n] = K[j,k]
            K_Freepres[n,m] = K[k,j]
            
    for i in range(np.size(FreeDof)):
            for j in range(np.size(FreeDof)):
                roundoff_error = 1.0E-12
                if i != j:
                    if (K_compute[i, j] - K_compute[j, i]) > roundoff_error:                    
                        print ('Check symmetric matrix')
                    else: continue
    
    u_compute = np.matmul((np.linalg.inv(K_compute)), \
    (p_compute + peq_compute -  np.matmul(K_Freepres, u_Prescribed)))    
    
    p_prescribed = np.matmul(K_presFree,u_compute) + \
    np.matmul(K_prescribed,u_Prescribed) - peq_prescribed
    
    for i, j in enumerate(FreeDof):
        u[j] = u_compute[i]
    for i,j in enumerate(idxPrescribedDof):
        p[j] = p_prescribed[i]
    
    return u,p,p_prescribed

u,p,p_prescribed = Solver(K, u, p, peq, bc)

print('Global displacements is U = \n', u)
print()
print('Global forces is P = \n', p)

