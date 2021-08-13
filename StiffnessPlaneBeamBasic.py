from input2d import *
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
