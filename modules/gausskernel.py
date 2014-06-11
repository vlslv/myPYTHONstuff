#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt 


def showModuleName():
    print __doc__

def gauss1d(media,sigma):
    x = np.linspace(media-5.*sigma,media+5.*sigma,60)
    y = (1./sigma/np.sqrt(2.*np.pi))*np.exp(-(x-media)**2/2./sigma**2)
    plt.plot(x,y)
    plt.grid()
    plt.show()
    return  
    
    
