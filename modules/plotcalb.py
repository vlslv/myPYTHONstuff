#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import dates

import numpy as np
import pylab
import pandas as pd
import random
from datetime import datetime, timedelta

from scipy.io.idl import readsav
from IPython.core.display import Image
import io
import base64
from IPython.display import HTML

import sunpy
from astropy.io import fits
from sunpy.sun import constants as con
from sunpy.net.helioviewer import HelioviewerClient
from sunpy.time import *
from sunpy.net import vso
from sunpy import lightcurve as lc
from sunpy.time import TimeRange
from sunpy.net import hek

def t_plot(fds,bi):
    """
    plota todos los valores 'adc' de la estrutura 'bi'
    contra el vector-tiempo contenido en 'fds'
    """
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    b1,b2,b3,b4,b5,b6 = ax.plot(fds,bi.adc,'-')
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_xlim(dates.datestr2num('2011/11/03 19:55:00'), dates.datestr2num('2011/11/03 20:00:00'))
    #ax.xaxis.set_major_locator(dates.MinuteLocator(30))
    #ax.set_ylim(bottom = 2.0e4,top=2.048e4)
    #ax.xaxis.set_major_locator(loc)
    plt.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5","beam 6"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    plt.show()
    return

def bi_plot(fds,bi):
    """
    plota todos los valores 'adc' de la estrutura 'bi'
    contra el vector-tiempo contenido en 'fds'
    """
    import matplotlib.pyplot as plt
    from matplotlib import dates
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(fds,bi.adc,'-')
    #ax.set_xlim(dates.datestr2num('2011/11/03 19:55:00'), dates.datestr2num('2011/11/03 20:00:00'))
    #ax.xaxis.set_major_locator(dates.MinuteLocator(30))
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_ylim(bottom = 2.0e4,top=2.048e4)
    #ax.xaxis.set_major_locator(loc)
    ax.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5" ,"beam 6"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    plt.show()
    return

def rs_plot(fds,bi):
    """
    plota todos los 'adcval' de la estrutura 
    'rs' contra el vector-tiempo contenido en fds
    """
    import matplotlib.pyplot as plt
    from matplotlib import dates
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(fds,bi.adcval,'-')
    #ax.set_xlim(dates.datestr2num('2011/11/03 19:55:00'), dates.datestr2num('2011/11/03 20:00:00'))
    #ax.xaxis.set_major_locator(dates.MinuteLocator(30))
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_ylim(bottom = 2.0e4,top=2.048e4)
    #ax.xaxis.set_major_locator(loc)
    legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5" ,"beam 6"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    plt.show()
    return

def caja(taxdict,ybot,ytop):
    """
    Boxplots del contenido del dictionario (taxdict)
    devuelto por la funcion 'taxonomia'. Se tiene que
    introducir los limites del eje vertical: 'ybot', 'ytop'.
    """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11.69*1.5,8.27*1.5))
    ax = fig.add_subplot(191)
    ax.set_title('cold/hot')
    ax.boxplot(taxdict['colval'],0,'b+')#ax.scatter(coldfds,coldadc[:,0],color='blue',marker='.')
    ax.boxplot(taxdict['hotval'],0,'r+')#ax.scatter(hotfds,hotadc[:,0],color='red',marker='.')
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(192)
    ax.set_title('trk')
    ax.boxplot(taxdict['trkval'],0,'g+')#ax.scatter(taufds,tauadc[:,0],color='skyblue')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(193)
    ax.set_title('map')
    ax.boxplot(taxdict['mapval'],0,'c+')#ax.scatter(taufds,tauadc[:,0],color='skyblue')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(194)
    ax.set_title('mpi')
    ax.boxplot(taxdict['mpival'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(195)
    ax.set_title('scn')
    ax.boxplot(taxdict['scnval'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(196)
    ax.set_title('sci')
    ax.boxplot(taxdict['scival'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(197)
    ax.set_title('tau')
    ax.boxplot(taxdict['tauval'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(198)
    ax.set_title('stl')
    ax.boxplot(taxdict['stlval'] ,1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(199)
    ax.set_title('unk')
    ax.boxplot(taxdict['unkval'],1,'y*')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    plt.show()
    return

def tx_plot(taxdict,str_key):
    """
    Plot de los adcval de cualquiera de las variables segregadas
    segun el diccionario (taxdict) devuelto por la funcion 'taxonomia'.
    Se tiene que entrar con la clave (str_key) para escoger el modo, la
    misma que aparece en los titulos del grafico porducido por 'caja'.
    """
    strfds = str_key+'time'
    stradc = str_key+'val'
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(taxdict[strfds],taxdict[stradc] ,'.-')
    legend([b1,b2,b3,b4,b5,b6], ["1","2","3","4","5","6"], loc='best')
    #ax.set_xlim(dates.datestr2num('2009/09/26 14:29:19'), dates.datestr2num('2009/09/26 14:29:50'))
    #ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_major_formatter(hfmt)
    #ax.set_ylim(bottom = 0)
    #ax.xaxis.set_major_locator(loc)
    plt.xticks(rotation=30)
    #plt.subplots_adjust(bottom=.3)
    #plt.legend(loc='best')
    plt.grid()
    plt.show()
    return

def st_plot(statdict):
    """
    Plot de la relacion lineal entre las estadisticas
    de cada canal de 212 GHz en relacion al canal 1 y
    del canal 5 (405 GHz) en relacion al canal 6
    """
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax1.plot(statdict['stat01'],statdict['stat01'],'bo')
    ax1.plot(statdict['stat01'],statdict['stat02'],'go')
    ax1.plot(statdict['stat01'],statdict['stat03'],'ro')
    ax1.plot(statdict['stat01'],statdict['stat04'],'co')
    legend(["1vs1","1vs2","1vs3","1vs4"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    
    ax2 = fig.add_subplot(1,2,2)
    ax2.plot(statdict['stat06'],statdict['stat06'],'mo')
    ax2.plot(statdict['stat06'],statdict['stat05'],'yo')
    legend(["6vs6","6vs5"], loc='best')
    plt.xticks(rotation=30)
    plt.grid()
    plt.show()
    return
