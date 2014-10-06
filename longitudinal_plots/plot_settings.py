'''
**Class to customize plot layout **

:Authors: **Helga Timko**

'''

import matplotlib.pyplot as plt


class PlotSettings(object):
    
    def __init__(self, linewidth=2, markersize=6, labelsize=18, 
                 fontfamily='sans-serif', fontweight='normal', dpi=100):
        '''
        Initialize custom plot formatting. For more options, see
        
        http://matplotlib.org/1.3.1/users/customizing.html
        '''

        self.lwidth = linewidth
        self.msize = markersize
        self.lsize = labelsize
        self.ffamily = fontfamily
        self.fweight = fontweight
        self.dpi = dpi
        
        # Ticksize
        self.tsize = self.lsize - 2
        
        
    def set_plot_format(self):
        
        # Set size of x- and y-grid numbers
        plt.rc('xtick', labelsize=self.tsize) 
        plt.rc('ytick', labelsize=self.tsize)
        
        # Set x- and y-grid labelsize and weight
        plt.rc('axes', labelsize=self.lsize)
        plt.rc('axes', labelweight=self.fweight)

        # Set linewidth for continuous, markersize for discrete plotting
        plt.rc('lines', linewidth=self.lwidth, markersize=self.msize)
        
        # Set figure resolution, font
        plt.rc('figure', dpi=self.dpi)  
        plt.rc('savefig', dpi=self.dpi)  
        plt.rc('font', family=self.ffamily)  
         
