# -*- coding: utf-8 -*-
"""

Created on Fri Nov 30 19:26:26 2016

@author: Jingying Hu
"""

import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import matplotlib.ticker as ticker

        
def plague(nGen = None, N = 214):
    
    mask = np.ones((10,10))
    mask[5,5] = 0
    """
    mask = np.array([[1, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 0, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1]])
    """
    
    # distribute residents outside=============================================
    class resident:
        def __init__(self, P_infect, P_recover, row, column):
            self.P_infect = np.random.choice((0, 1))
            self.P_recover = np.random.choice((0, 1))
            self.row = row
            self.column = column 
            
    outside_residents = []
    for i in range(14):
        outside_residents.append(resident(P_infect = np.random.choice((0, 1)),
                                          P_recover = np.random.choice((0, 1)),
                                          row = 0, column = 0))
    for resident in outside_residents:
        resident.row = np.random.choice(range(100))
        if 15 < resident.row < 84:
            list1 = np.arange(15)
            list2 = np.arange(84,100)
            list0 = np.concatenate((list1,list2))
            resident.column = np.random.choice(list0)
        else:
            resident.column = np.random.choice(range(100))
    z1 = np.zeros((14,2))
    for i in range(14):
        row = outside_residents[i].row
        column = outside_residents[i].column
        z1[i,0] = column
        z1[i,1] = row
        
    # distribute resides in the middle=========================================
    class resident:    
        def __init__(self, P_infect, P_recover, row, column):
            self.P_infect = np.random.choice((0, 1))
            self.P_recover = np.random.choice((0, 1))
            self.row = row
            self.column = column     
    
    middle_residents = []
    for i in range(85):
        middle_residents.append(resident(P_infect = np.random.choice((0, 1)),
                                          P_recover = np.random.choice((0, 1)),
                                          row = 0, column = 0))
    for resident in middle_residents:
        resident.row = np.random.choice(range(15, 84))
        if 37 < resident.row < 63:
            list1 = np.arange(15,37)
            list2 = np.arange(63,84)
            list0 = np.concatenate((list1,list2))
            resident.column = np.random.choice(list0)
        else:
            resident.column = np.random.choice(range(15,84))
            
    z2 = np.zeros((85,2))        
    for i in range(85):
        row = middle_residents[i].row
        column = middle_residents[i].column
        z2[i,0] = column
        z2[i,1] = row
                           
    # distribute resides in the center========================================= 
    class resident:    
        def __init__(self, P_infect, P_recover, row, column):
            self.P_infect = np.random.choice((0, 1))
            self.P_recover = np.random.choice((0, 1))
            self.row = row
            self.column = column        
            
    center_residents = []
    for i in range(115):
         center_residents.append(resident(P_infect = np.random.choice((0, 1)),
                                          P_recover = np.random.choice((0, 1)),
                                          row = 0, column = 0))
    for resident in center_residents:    
        resident.row = np.random.choice(range(37,63))
        resident.column = np.random.choice(range(37,63))
        
    z3 = np.zeros((115,2))
    for i in range(115):
        row = center_residents[i].row
        column = center_residents[i].column
        z3[i,0] = column
        z3[i,1] = row

    #concatenate z1,z2,z3
    coordinates = np.concatenate((z1,z2,z3))  # coordinates of each resident      
    # find sick
    d = np.random.choice(range(N))
    sick = coordinates[d,:] # coordinate of the sick resident
    
    # framework
    n = 100
    z = np.zeros((n, n), dtype = int)
   
    # assign value '1' to the healthy
    for i in range(N):
        a = int(coordinates[i][1])
        b = int(coordinates[i][0])
        z[a, b] = 2 # healthy, green
        
    # assign value '2' to the sick     
    z[int(sick[1]),int(sick[0])] = 4   #sick, red  
    
    # colormap      
    cdict = {'red':   ((0.00, 1.00, 1.00),
                       (0.33, 0.00, 0.00),
                       (0.67, 1.00, 1.00),
                       (1.00, 0.47, 0.47)),
             'green': ((0.00, 1.00, 1.00),
                       (0.33, 0.80, 0.80),
                       (0.67, 0.40, 0.40),
                       (1.00, 0.47, 0.47)),
             'blue':  ((0.00, 1.00, 1.00),
                       (0.33, 0.00, 0.00),
                       (0.67, 0.00, 0.00),
                       (1.00, 0.47, 0.47))}
    colormap = colors.LinearSegmentedColormap('mycolors', cdict, 256)
    fig, ax = plt.subplots()
    plt.axis('scaled')
    plt.axis([0, n, 0, n]) 
    cplot = plt.pcolormesh(z, cmap = colormap, vmin = 0, vmax = 6)
    plt.title('plague spreading')
    
    # animation updata function
    def update(i):  
        """
        States: empty (0) 'white', the healthy (2) 'green', the sick (4) 'red', 
        the dead (6) 'grey', the immuned (1) 'light green'
        Rules:
          rule 1: green tree with one or more burning neighbors bursts into flame
          rule 2: burning tree becomes a burnt stump
        """
        n = ndimage.generic_filter(z == 4, np.sum, footprint = mask, mode = 'constant', output = int)
        r1 = (z == 4) & (n > 0)
        r2 = (z == 4) & (np.random.binomial(1,0.5) == 1) 
        z[r1] = 4 # infected
       # r3 = (z == 4) & (?)
        z[r2] = 6 # dead, grey
       # z[r3] = 1 # immuned light green
        cplot.set_array(z.ravel()) #set_array requires a 1D array
        plt.title('plague spreading')
        return cplot
    
    def genner():
        i = 1
        while True:
            if (z == 2).any(): #if any trees are burning        
                yield i
                i += 1
            else:
                return    
    
    if nGen is None:    
        nGen = genner()
    
    #Time to put those midi-chlorians to work, Ani    
    anakin = animation.FuncAnimation(fig, update,
                           frames = nGen, 
                           fargs = (),
                           #init_func = initialize,
                           interval = 400, 
                           blit = False, #cplot is not iterable, so blit must be False
                           repeat = False)
    return anakin

        
  
        


