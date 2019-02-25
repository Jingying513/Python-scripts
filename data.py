# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 23:09:33 2019

@author: jy940
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 19:10:57 2018

@author: Administrator
"""
"""
0: Acidification
1: Ecotoxicity
2: Eutrophication
3: Global Warming
4: carcinogenic
5: non-carcinogenic
6: Ozone Depletion
7: Photochemical ozone formation
8: Resource Depletion
9: Respiratory effects    
"""
# grab information from openLCA Simulation result
import pandas as pd
import os
import numpy as np

#dfs = pd.read_excel('C:/Users/jy940/iCloudDrive/MOYYEF.xlsx', index_col='Impact category', sheet_name='Impact Assessment')
dfs = []
c
for root, dirs, files in os.walk('D:/RECIPE (H) result'):
    print(files)
    metald = np.zeros((len(files),1)) 
    for i in range(len(files)):
        filename = 'D:/RECIPE (H) result/' + files[i] 
        dfs = pd.read_excel(filename,
                      index_col = 'Impact category',
                      sheet_name = 'Impacts')
        metald[i] = dfs.iloc[9, 3]
        
np.savetxt('C:/Users/jy940/Desktop/Metal_depletion.txt', metald)

