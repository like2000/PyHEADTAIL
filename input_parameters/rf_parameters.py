'''
Created on 19 juin 2014

Module gathering and processing all the RF parameters to be given to the other modules

@author: alasheen
'''

class RF_parameters(object):
    '''
    Object containing all the RF parameters
    '''

    def __init__(self, n_sections):

        self.n_sections = n_sections #: Number of sections of the ring