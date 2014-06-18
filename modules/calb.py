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

def hola():
    print 'Hello World, this my calibration module'
    return

def lee_sst(file_place):
    """
    lee archivos '/place/of/file/bi1yymmdd.save' o 
    '/place/of/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (biors).
    """
    pathplusfile = file_place 
    biors = readsav(pathplusfile,python_dict=False)
    yyyy = np.int(pathplusfile[37:39])+2000
    mm = np.int(pathplusfile[39:41])
    dd = np.int(pathplusfile[41:43])
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    thatday = obsday + biors.time/1.e4
    dts = map(datetime.fromtimestamp,thatday)
    fds = dates.date2num(dts)
    return fds,biors

def modos(biors):
    """
    Devuelve en pantalla la lista de modos existentes
    en el archivo 'biors' y los guarda en el diccionario
    'opmodict'.
    """
    modos=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,40,50,55,99]
    opmodict={'0':0,'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0,'9':0,'10':0,'11':0,\
          '12':0,'13':0,'14':0,'15':0,'16':0,'17':0,'18':0,'19':0,'20':0,'21':0,'22':0,\
          '23':0,'40':0,'50':0,'55':0,'99':0}
    for i in modos:
        xx=np.where(biors.opmode == i)
        tama = np.size(xx)
        if (tama >  0):
            opmodict[str(i)] = tama
    for i in modos:
        xx=np.where(biors.opmode == i)
        tama = np.size(xx)
        if (tama >  0):
            print i, tama
    return opmodict

def taxonomia(fdsh,bh):
    """
    Separa los datos adcval y los correspondientes 
    vectores-tiempo segun los modos mas comunes del
    'opmode'. Guarda los datos segregados en el dic-
    cionario 'taxdict'. Funciona solo con archivos de
    tipo rs1yymmdd.hh00.save
    """
    colfds = fdsh[np.where(bh.target/32 == 1)]
    coladc = bh.adcval[np.where(bh.target/32 == 1)]
    hotfds = fdsh[np.where(bh.target/32 == 2)]
    hotadc = bh.adcval[np.where(bh.target/32 == 2)]
    trkadc = bh.adcval[np.where(bh.opmode == 0)]
    trkfds = fdsh[np.where(bh.opmode == 0)]
    mapadc = bh.adcval[np.where(bh.opmode == 2)]
    mapfds = fdsh[np.where(bh.opmode == 2)]
    mpiadc = bh.adcval[np.where(bh.opmode == 4)]
    mpifds = fdsh[np.where(bh.opmode == 4)]
    scnadc = bh.adcval[np.where(bh.opmode == 5)]
    scnfds = fdsh[np.where(bh.opmode == 5)]
    sciadc = bh.adcval[np.where(bh.opmode == 9)]
    scifds = fdsh[np.where(bh.opmode == 9)]
    tauadc = bh.adcval[np.where(bh.opmode == 10)]
    taufds = fdsh[np.where(bh.opmode == 10)]
    stladc = bh.adcval[np.where(bh.opmode == 50)]
    stlfds = fdsh[np.where(bh.opmode == 50)]
    unkadc = bh.adcval[np.where(bh.opmode == 99)]
    unkfds = fdsh[np.where(bh.opmode == 99)]
    taxdict={\
             'coltime':colfds,'colval':coladc\
             ,'hottime':hotfds,'hotval':hotadc\
             ,'trktime':trkfds,'trkval':trkadc\
             ,'maptime':mapfds,'mapval':mapadc\
             ,'mpitime':mpifds,'mpival':mpiadc\
             ,'scntime':scnfds,'scnval':scnadc\
             ,'scitime':scifds,'scival':sciadc\
             ,'tautime':taufds,'tauval':tauadc\
             ,'stltime':stlfds,'stlval':stladc\
             ,'unktime':unkfds,'unkval':unkadc\
             }
    return taxdict

def colecta(taxdict):
    """
    Escoge percentiles 25,50,75  de cada modo esco-
    gido por 'taxonomia' y compone un diccionario de
    estadisticas 'origdict'.
    """
    origdict = {}
    if len(taxdict['colval']) > 0:
        coldadc = taxdict['colval']
        statcol = [\
        np.percentile(coldadc,25,axis=0),np.percentile(coldadc,50,axis=0),np.percentile(coldadc,75,axis=0)]
        extra = {'colstat':statcol}
        origdict.update(extra)
    if len(taxdict['hotval']) > 0:
        hotadc = taxdict['hotval']
        stathot = [\
        np.percentile(hotadc,25,axis=0),np.percentile(hotadc,50,axis=0),np.percentile(hotadc,75,axis=0)]
        extra = {'hotstat':stathot}
        origdict.update(extra)
    if len(taxdict['trkval']) > 0:
        trkadc = taxdict['trkval']
        stattrk = [\
        np.percentile(trkadc,98,axis=0),np.percentile(trkadc,99,axis=0),np.percentile(trkadc,100,axis=0)]
        extra = {'trkstat':stattrk}
        origdict.update(extra)
    if len(taxdict['mapval']) > 0:
        mapadc = taxdict['mapval']
        statmap = [\
        np.percentile(mapadc,50,axis=0),np.percentile(mapadc,50,axis=0),np.percentile(mapadc,75,axis=0)]
        extra = {'mapstat':statmap}
        origdict.update(extra)
    if len(taxdict['mpival']) > 0:
        mpiadc = taxdict['mpival']
        statmpi = [\
        np.percentile(mpiadc,50,axis=0),np.percentile(mpiadc,50,axis=0),np.percentile(mpiadc,75,axis=0)]
        extra = {'mpistat':statmpi}
        origdict.update(extra)
    if len(taxdict['scnval']) > 0:
        scnadc = taxdict['scnval']
        statscn = [\
        np.percentile(scnadc,25,axis=0),np.percentile(scnadc,50,axis=0),np.percentile(scnadc,75,axis=0)]
        extra = {'scnstat':statscn}
        origdict.update(extra)
    if len(taxdict['scival']) > 0:
        sciadc = taxdict['scival']
        statsci = [\
        np.percentile(sciadc,25,axis=0),np.percentile(sciadc,50,axis=0),np.percentile(sciadc,75,axis=0)]
        extra = {'scistat':statsci}
        origdict.update(extra)
    if len(taxdict['tauval']) > 0:
        tauadc = taxdict['tauval']
        stattau = [\
        np.percentile(tauadc,25,axis=0),np.percentile(tauadc,50,axis=0),np.percentile(tauadc,75,axis=0)]
        extra = {'taustat':stattau}
        origdict.update(extra)
    if len(taxdict['stlval']) > 0:
        stladc = taxdict['stlval']
        statstl = [\
        np.percentile(stladc,25,axis=0),np.percentile(stladc,50,axis=0),np.percentile(stladc,75,axis=0)]
        extra = {'stlstat':statstl}
        origdict.update(extra)
    if len(taxdict['unkval']) > 0:
        unkadc = taxdict['unkval']
        statunk = [\
        np.percentile(unkadc,25,axis=0),np.percentile(unkadc,50,axis=0),np.percentile(unkadc,75,axis=0)]
        extra = {'unkstat':statunk}
        origdict.update(extra)
    return origdict

def pseudocal(origidict):
    """
    Se trata de colocar las estadisticas de obtenidas de
    'colecta' y almacenados en el dicccionario 'origdict'
    para hacer la pseudo-calibracion que iguala canales.
    """
    statdict = {}
    stat1 = [0.0]
    stat2 = [0.0]
    stat3 = [0.0]
    stat4 = [0.0]
    stat5 = [0.0]
    stat6 = [0.0]
    for key in origidict.keys():
        var1 = [origidict[key][0][0],origidict[key][1][0],origidict[key][2][0]]
        stat1.extend(var1)
        var2 = [origidict[key][0][1],origidict[key][1][1],origidict[key][2][1]]
        stat2.extend(var2)
        var3 = [origidict[key][0][2],origidict[key][1][2],origidict[key][2][2]]
        stat3.extend(var3)
        var4 = [origidict[key][0][3],origidict[key][1][3],origidict[key][2][3]]
        stat4.extend(var4)
        var5 = [origidict[key][0][4],origidict[key][1][4],origidict[key][2][4]]
        stat5.extend(var5)
        var6 = [origidict[key][0][5],origidict[key][1][5],origidict[key][2][5]]
        stat6.extend(var6)
    del stat1[0],stat2[0],stat3[0],stat4[0],stat5[0],stat6[0]    
    statdict = {'stat01':stat1,'stat02':stat2,'stat03':stat3,'stat04':stat4,'stat05':stat5,'stat06':stat6}
    return statdict 
