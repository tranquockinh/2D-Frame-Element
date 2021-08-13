import numpy as np
## ------------------------------------------- Please input data -------------------------------------------------
# Defind the structure coordinates
ndcoords = 12 * np.array([[0.0, 0.0],
                           [0.0, 10],
                           [10, 10],
                           [10, 0.0]])
# Defining nodal connection
elems = np.array([[0, 1],
                   [1, 2],
                   [2, 3]])                   
# Boundary conditions
bc = np.array([[0, 0, 0],
                [0, 1, 0],
                [0, 2, 0],
                [3, 0, 0],
                [3, 1, 0],
                [3, 2, 0]], dtype=int)  
# Material properties
'''
## Use this in new problem
b = np.array([2.5, 2.5, 2.5]) ## width
h = np.array([4, 4, 4]) ## height
A = b * h ## area
ymax = h / 2 ## length from center to outermost in y-direction
E = np.array([30E+06, 30E+06, 30E+06]) ## elatic modulus
I = b * h**3 / 12
## This can be a loop to automatically append the data into eparray matrix
eparray = np.matrix([[E[0]*A[0], E[0]*I[0]],
                     [E[1]*A[1], E[1]*I[1]],
                     [E[2]*A[2], E[2]*I[2]]])
'''
## The example below is taken from the book: "The first course in FEM, DARYL. L. Logan, p.239-244"
#                        EA       EI
E = 30E6
eparray = np.matrix([[30E6*10, 30E6*200],
                     [30E6*10, 30E6*100],
                     [30E6*10, 30E6*200]])
# Concentrated load
cloads = np.matrix([[1,     10000,      0,      0],
                    [2,     0,          0,      5000]])
# Distributed load
qloads = np.matrix([[1,      0,      0,      0],
                    [2,      0,      0,      0]])        
## ------------------------------------------- Ending input data -------------------------------------------------
        
##=========PLEASE CHANGE THE INITIALIZATION BELOW ONLY WHEN NEEDED
# Initialization
ndim = 3
numnodes = np.size(ndcoords, 0) # number of elements in the first column
ndof = ndim * numnodes
nelems = np.size(elems, 0)
K = np.zeros((ndof, ndof),dtype=float)
ke = np.zeros((ndim * 2, ndim * 2),dtype=float)
##print ('K is \n {} \n and \n ke is \n {}'.format(K, ke)) # Checking
p = np.zeros((ndof, 1)) # to ensure the force vector equals ndof
peq = np.zeros((ndof, 1))
pe = np.zeros((ndim * 2, 1))
u = np.zeros((ndof,1))
## suppose A = 10 (inxin) square cross-section
ymax = np.sqrt(10) / 2


