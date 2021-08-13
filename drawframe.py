import numpy as np
import matplotlib.pyplot as plt


from input2d import ndcoords, elems, bc
'''
def input(ndcoords, elems, bc):
    coords  = ndcoords
    connect = elems
    bcs = bc
    return coords, connect, bcs
coords, connect, bcs = input(ndcoords, elems, bc)
print(connect)
'''
'''
ndcoords = np.array([[0.0, 0.0],
                   [4.0, 0.0],
                   [8.0, 0.0],
                   [12.0, 0.0],
                   [4.0, -6.0],
                   [8.0, -8.0],
                   [8.0, -4.0]])

elems = np.array([[0, 1],
                    [1, 2],
                    [2, 3],
                    [4, 1],
                    [5, 6],
                    [6, 2]])
                               
bc = np.array([[0, 2, 0],
                [3, 1, 0],
                [4, 0, 0],
                [4, 1, 0],
                [4, 2, 0],
                [5, 0, 0],
                [5, 1, 0],
                [5, 2, 0]], dtype=int)

'''
# Not to change the code below, ca be changed if need to modify the code
def draw_frame(ndcoords, elems, bc):
    prescribed = np.unique(bc[:,0])    
    component = np.empty((len(prescribed),1),dtype=object)
                
    bc_trans = np.zeros((len(prescribed),3),dtype='int')
    bc_trans = np.concatenate((prescribed[:, np.newaxis],bc_trans),axis=1)  


    for i,j in enumerate(prescribed):
        component[i,0]=np.array(np.where(bc[:,0]==j)).flatten()
        
        if len(component[i,0])==1:
            node_number = component[i,0][0]
            #locate = np.array(np.where[bc[:,0]==node_number]).flatten()
            if bc[node_number,1]==0:
                bc_trans[i,1:] = np.array([0,1,1])
            elif bc[node_number,1]==1:
                bc_trans[i,1:] = np.array([1,0,1])
            else: 
                bc_trans[i,1:] = np.array([1,1,0])
        elif len(component[i,0])==2:
            first = component[i,0][0]
            second = component[i,0][1]
            if (bc[first,1]+bc[second,1])==1:
                bc_trans[i,1:] = np.array([0,0,1])
            elif (bc[first,1]+bc[second,1])==2:
                bc_trans[i,1:] = np.array([0,1,0])
            else:
                bc_trans[i,1:] = np.array([1,0,0])
        else:
            bc_trans[i,1:] = np.array([0,0,0])

    print(bc_trans)
    ##===============================================
    ## Note: '0' means 'to be prescribed'
    ##       '1' means free
    ## Example
    '''
    ## The first column is the node numbers to be prescribed        
    bc = np.array([[0,1,0,1], ## x, and theta are free, y is constrained
                   [3,1,0,1], ## x, and theta are free, y is constrained
                   [4,0,0,0], ## x, and theta and y are constrained
                   [5,0,0,0]], dtype = 'int')  ## x, and theta and y are constrained
    '''
    ##===============================================

                   
    node_array = np.arange(0,len(ndcoords[:,0]),1)
    elem_array = np.arange(0,len(elems[:,0]),1)
    ndcoords = np.concatenate((node_array[:, np.newaxis],ndcoords),axis=1)
    elems = np.concatenate((elem_array[:, np.newaxis],elems),axis=1)               
    numnode = len(ndcoords[:,0]) ## number of nodes
    numel = len(elems[:,0]) ## number of elements
    numcon = 2 ## number of nodes that are limited from translations or rotation
             
    fig, ax = plt.subplots()
    start = np.zeros((np.size(elems,0),1),dtype=object)
    end = np.zeros((np.size(elems,0),1),dtype=object)
    Nodelabels = ndcoords[:,0]
    Elelabels = elems[:,0]

    for i in range(numel):
        
        ## Draw frame
        left = elems[i,1]
        right = elems[i,2]
        locateleft = np.array(np.where(ndcoords[:,0]==left)).flatten()
        locateright = np.array(np.where(ndcoords[:,0]==right)).flatten()
        start[i,0] = ndcoords[locateleft[0],1:]
        end[i,0] = ndcoords[locateright[0],1:]
        ax.plot([start[i,0][0],end[i,0][0]],[start[i,0][1],end[i,0][1]],'ks-',linewidth=3)
        
        ## Assign node numbers
        for j, addNodeNum in enumerate(Nodelabels):
            ax.annotate(int(addNodeNum), (ndcoords[j,1],ndcoords[j,2]),\
            textcoords='offset points',xytext=(8,8),ha='center',fontsize=12,color='blue')
            fig.patch.set_visible('False')
        
        ## Assign element numbers
        for j, addElemNum in enumerate(Elelabels):
            connection = elems[j,1:]
            leftpoint = ndcoords[np.where(ndcoords[:,0]==connection[0]),0]
            rightpoint = ndcoords[np.where(ndcoords[:,0]==connection[1]),0]
            placeX = ndcoords[int(leftpoint),1] + 1/2 * (ndcoords[int(rightpoint),1]-ndcoords[int(leftpoint),1])
            placeY = ndcoords[int(leftpoint),2] + 1/2 * (ndcoords[int(rightpoint),2]-ndcoords[int(leftpoint),2])
           
            ax.annotate(addElemNum, (placeX,placeY), \
            textcoords='offset points',xytext=(8,8),ha='center',fontsize=12,\
            bbox=dict(facecolor='none', edgecolor='blue', pad=2.5))
        
        ## Draw constrains
        for node,j in enumerate(bc_trans[:,0]):
            if np.sum(bc_trans[node,1:])==0:
                xcoords = ndcoords[j,1]
                ycoords = ndcoords[j,2]
                ax.plot(xcoords,ycoords,'rs',markersize=12, fillstyle='none')    
            elif np.sum(bc_trans[node,1:])==1: 
                if ((bc_trans[node,1])!=0 and (bc_trans[node,2] and bc_trans[node,3])==0):
                        xcoords = ndcoords[j,1]
                        ycoords = ndcoords[j,2]
                        ax.plot(xcoords,ycoords,'r_',markersize=20)
                        ax.plot(xcoords,ycoords,'ro',markersize=10, fillstyle='none')
                elif ((bc_trans[node,2])!=0 and (bc_trans[node,1] and bc_trans[node,3])==0):
                    xcoords = ndcoords[j,1]
                    ycoords = ndcoords[j,2]
                    ax.plot(xcoords,ycoords,'r|',markersize=20)
                else:
                    xcoords = ndcoords[j,1]
                    ycoords = ndcoords[j,2]
                    ax.plot(xcoords,ycoords,'r^',markersize=15, fillstyle='none')
                    ax.plot(xcoords,ycoords,'r.',markersize=10)
            else:  
                if ((bc_trans[node,1] and bc_trans[node,3])!=0 and (bc_trans[node,2])==0):
                        xcoords = ndcoords[j,1]
                        ycoords = ndcoords[j,2]
                        ax.plot(xcoords,ycoords,'ro',markersize=10, fillstyle='none')
                        ax.plot(xcoords,ycoords,'r.',markersize=10)
                        ax.plot(xcoords,ycoords,'r_',markersize=20)
                elif ((bc_trans[node,1] and bc_trans[node,2])!=0 and (bc_trans[node,3])==0):
                        xcoords = ndcoords[j,1]
                        ycoords = ndcoords[j,2]
                        #ax.plot(xcoords,ycoords,'ro',markersize=10, fillstyle='none')
                        ax.plot(xcoords,ycoords,'r+',markersize=20)
                        ax.plot(xcoords,ycoords,'ro',markersize=10, fillstyle='none')
                else:
                    xcoords = ndcoords[j,1]
                    ycoords = ndcoords[j,2]
                    ax.plot(xcoords,ycoords,'ro',markersize=10, fillstyle='none')
                    ax.plot(xcoords,ycoords,'r.',markersize=10)
                    ax.plot(xcoords,ycoords,'r|',markersize=20)

    ax.axis('off')  
    ax.set_title('Frame structure',loc='left')
    plt.show()   
    
    return ax

ax = draw_frame(ndcoords, elems, bc)  