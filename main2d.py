from input2d import *
from drawframe import *
from StiffnessPlaneBeamBasic import StiffnessPlaneBeamBasic
from StiffnessPlaneBeam import StiffnessPlaneBeam
from AssemblePlaneBeam import AssemblePlaneBeam
from AssemblePlaneBeamForces import AssemblePlaneBeamForces
from Solver import Solver

p = AssemblePlaneBeamForces(p, cloads)

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
        ## Equivalent loads
        peq = AssemblePlaneBeamForces(peq, eqvcloads)
    else:        
        ## Call the function of stiffness matrix
        Ke,T, kle, L, b = StiffnessPlaneBeamBasic(ex, ey, ep) 
        
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
    Global element stiffness matrix element {}=\n{}\nT{}=\n {}'.format(ThisElement,Ke,ThisElement,T)) 
    print()
    
u,p,u_reshaped,p_reshaped = Solver(K, u, p, peq, bc)

print('---------------------------OUTPUTS---------------------------')
print('Global displacements is U = \n', u)
print()
print('Global forces is P = \n', p)
# Compute element forces, stresses
u_global = np.zeros((ndim*2, len(elems[:,0])))
u_local = np.zeros((ndim*2, len(elems[:,0])))
p_global = np.zeros((ndim*2, len(elems[:,0])))
p_local = np.zeros((ndim*2, len(elems[:,0])))
bending_stress = np.zeros((2, nelems))
axial_stress = np.zeros((nelems))
## Element displacement and local element forces
for j in range(nelems):
    LeftNode = elems[j, 0]
    RightNode = elems[j, 1]
    ex = [ndcoords[LeftNode, 0], ndcoords[RightNode, 0]]
    ey = [ndcoords[LeftNode, 1], ndcoords[RightNode, 1]]
    ep = eparray[j, :]
    found = 0 
    for i in range(np.size(qloads,0)):
        if j == qloads[i,0]:
            found = 1
            break
    if found:
        Ke, T, kle, pe, L, b = StiffnessPlaneBeam(ex, ey, ep, qloads[i, 1:])
        pe = np.array(pe.T).flatten() ## Fattenning pe to an array
        peLeft = np.append([LeftNode], pe[:3])
        peRight = np.append([RightNode], pe[3:])
        eqvcloads = np.matrix([peLeft, peRight])
        peq = AssemblePlaneBeamForces(peq, eqvcloads)
    else:        
        Ke,T, kle, L, b = StiffnessPlaneBeamBasic(ex, ey, ep) 
        
    u_global[:,j] = np.concatenate((u_reshaped[:,LeftNode],u_reshaped[:,RightNode]), axis=0)
    p_global[:,j] = np.concatenate((p_reshaped[:,LeftNode],p_reshaped[:,RightNode]), axis=0)
    u_local[:,j][:, np.newaxis] = T * u_global[:,j][:, np.newaxis]
    p_local[:,j][:, np.newaxis] = kle * T * u_global[:,j][:, np.newaxis]
    
    ## Compute stress and moment of each element
    L = np.sqrt((ex[1]-ex[0])**2+(ey[1]-ey[0])**2)
    
    stress_coeff_matrix = np.array([[6/L,-2/(L)],
                                    [6/L, 2/(L)]])
                                   
    disp_vector = np.array([[u_global[4,j] - u_global[1,j]],
                            [2*u_global[2,j] + u_global[5,j]]])
                            
    bending_stress[0,j] = ((E*ymax) * np.matmul(stress_coeff_matrix, disp_vector))[0]
    bending_stress[1,j] = ((E*ymax) * np.matmul(stress_coeff_matrix, disp_vector))[1]
    
    ## Compute eleemt axial stress
    axial_stress[j] = E * np.matmul(np.array([-1/L, 1/L]), np.array([[u_global[0,j]],[u_global[3,j]]]))

print()
print('Displacement at each node:\n',u_local)
print()
print('Element forces of each member:\n',p_local)
print()
print('Element bending stresses: with (+/-) sign\n', bending_stress)
print()
print('Element axial stresse:\n', axial_stress)
