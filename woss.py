# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 18:26:02 2014

@author: pedro.correia
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def do_linear_transform(x,minimum,maximum):
    dmin = x.min()
    dmax = x.max()
    data = (x-dmin)*(maximum-minimum)/(dmax-dmin)+minimum
    return data

def interpolate(x,y):
    fx = interp1d(y,x,kind='cubic')
    xf = np.float_(fx(np.linspace(y.min(),y.max(),1000)))
    return xf

def do_pearson_correlation(a,b):
    return np.sum((a-a.mean())*(b-b.mean()))/np.sqrt(np.sum((a-a.mean())**2)*np.sum((b-b.mean())**2))
    
def do_reflective_similarity(a,b):
    return np.sum(a*b)/np.sqrt(np.sum(a**2)*np.sum(b**2))
    
def do_quasi_correlation(a,b):
    return 2*np.sum(a*b)/(np.sum(a**2)+np.sum(b**2))
    
def do_squared_quasi_correlation(a,b):
    return (2*np.sum(a*b)/(np.sum(a**2)+np.sum(b**2)))**2

def do_spearman_correlation(a,b):
    return np.sum((a-a.mean())*(b-b.mean()))/np.sqrt(np.sum((a-a.mean())**2)*np.sum((b-b.mean())**2))
    
def do_cosine_correlation(a,b):
    return np.sum(a*b)/np.sqrt(np.sum(a**2)*np.sum(b**2))

def convolve_trace(trace,wavelet):
    reflectivity = trace[:,-1].copy()
    for i in xrange(trace.shape[0]-1):
        reflectivity[i]=(reflectivity[i+1]-reflectivity[i])/(reflectivity[i+1]+reflectivity[i])
    reflectivity[-1]=0
    synthetic = trace[:,-1].copy()
    synthetic[:] = 0
    h_size=(wavelet.shape[0]-1)/2
    for i in xrange(trace.shape[0]):
        if i-h_size<0:
            wa=h_size-i
            a=0
        else:
            wa=0
            a=i-h_size
        if i+h_size>trace.shape[0]:
            wb=h_size+i-trace.shape[0]
            b=trace.shape[0]
        else:
            wb=2*h_size+1
            b=i+h_size
        synthetic[a:b]=synthetic[a:b]+reflectivity[i]*wavelet[wa:(2*h_size-wb)]
    return synthetic
    
def convolve_trace2(trace,wavelet):
    reflectivity = trace[:,-1].copy()
    for i in xrange(trace.shape[0]-1):
        reflectivity[i]=(reflectivity[i+1]-reflectivity[i])/(reflectivity[i+1]+reflectivity[i])
    reflectivity[-1]=0
    synthetic = trace[:,-1].copy()
    synthetic[:] = 0
    h_size=(wavelet.shape[0]-1)/2
    for i in xrange(trace.shape[0]):
        if i-h_size<0:
            wa=h_size-i
            a=0
        else:
            wa=0
            a=i-h_size
        if i+h_size>trace.shape[0]:
            wb=h_size+i-trace.shape[0]
            b=trace.shape[0]
        else:
            wb=2*h_size+1
            b=i+h_size
        synthetic[a:b]=synthetic[a:b]+reflectivity[i]*wavelet[wa:(2*h_size-wb),-1]
    return synthetic

def do_ricker_wavelet(f, size, dt=1):
    t = np.int_(np.linspace(-size, size, (2*size+1)/dt))
    y = (1.0 - 2.0*(np.pi**2)*(f**2)*(t**2)) * np.exp(-(np.pi**2)*(f**2)*(t**2))
    data = np.hstack((t[:,np.newaxis],y[:,np.newaxis]))
    return data

def woss(well_paths,wavelet_path,wavelet_size,seismic_path,ites=1000,f=0.3,size=29):
    well_book = {}
    seismic_book = {}
    seismic = np.load(seismic_path)
    c=0
    initial_wavelet = np.loadtxt(wavelet_path)
    ricker_wavelet = do_ricker_wavelet(f,size)
    awavelet = np.zeros((ricker_wavelet.shape[0],2))
    awavelet[:,0] = np.linspace(initial_wavelet[:,0].min(),initial_wavelet[:,0].max(),2*size+1)
    awavelet[:,1] = do_linear_transform(ricker_wavelet[:,1],initial_wavelet[:,1].min(),initial_wavelet[:,1].max())
    initial_wavelet = awavelet.copy()
    
    final_wavelet = np.zeros(initial_wavelet.shape)
    final_wavelet[:,-1] = initial_wavelet[:,-1]
    for i in well_paths:
        well_book[c] = np.loadtxt(i)
        seismic_book[c] = seismic[well_book[c][:,0].astype('int32'),well_book[c][:,1].astype('int32'),well_book[c][:,2].astype('int32')]
        c=c+1
    cc0 = 0
    for t in well_book.keys():
        tt = convolve_trace2(well_book[t],final_wavelet)
        cc0 = cc0 + do_quasi_correlation(tt,seismic_book[t])
    ccs = np.float(cc0)
    ccf = np.float(cc0)
    nt = []
    for i in xrange(ites):
        c = np.random.randint(0,initial_wavelet.shape[0])
        p = np.random.triangular(initial_wavelet[:,-1].min(),final_wavelet[c,-1],initial_wavelet[:,-1].max())
        appex = final_wavelet[c,-1]
        final_wavelet[c,-1]=p
        cc = 0
        for t in well_book.keys():
            tt = convolve_trace2(well_book[t],final_wavelet)
            cc = cc + do_quasi_correlation(tt,seismic_book[t])
        if cc>cc0:
            cc0 = np.float(cc)
            ccf = np.float(cc0)
            nt.append(i)
        else:
            final_wavelet[c,-1] = appex 
        print 'START: ',ccs/len(well_book.keys()), ' - ONGOING: ',ccf/len(well_book.keys())
        print 'DONE ITERATION ',i
    print '#############################################'
    print 'INITIAL MEAN CORRELATION WAS: ',ccs/len(well_book.keys())
    print 'FINAL MEAN CORRELATION IS: ',ccf/len(well_book.keys())
    plt.plot(initial_wavelet[:,1],initial_wavelet[:,0],color='red',linewidth=3,label='initial wavelet')
    plt.plot(final_wavelet[:,1],initial_wavelet[:,0],color='black',linewidth=3,label='optimized wavelet')
    plt.plot(interpolate(final_wavelet[:,1],initial_wavelet[:,0]),np.linspace(initial_wavelet[:,0].min(),initial_wavelet[:,0].max(),1000),color='green',linewidth=3,label='interpolation')
    plt.fill_betweenx(initial_wavelet[:,0],initial_wavelet[:,1],final_wavelet[:,1],color='pink',alpha=0.3)    
    plt.xlim(initial_wavelet[:,1].min(),initial_wavelet[:,1].max())
    plt.ylim(initial_wavelet[:,0].min(),initial_wavelet[:,0].max())
    x=[]
    y=[]
    z=[]
    for i in xrange(initial_wavelet[:,0].shape[0]):
        y.append(initial_wavelet[i,0])
        y.append(initial_wavelet[i,0])
        z.append(initial_wavelet[i,1]-final_wavelet[i,1])
        z.append(initial_wavelet[i,1]-final_wavelet[i,1])
        x.append(initial_wavelet[i,1])
        x.append(final_wavelet[i,1])
    plt.legend()
    plt.grid()
    plt.show()
    
######### THIS IS HOW I WOULD LAUNCH THIS ALGORITHM ############
#well_paths = ['logs_AI_1.prn','logs_AI_2.prn']
#wavelet_path = 'wavelet.txt'
#wavelet_size = 19
#seismic_path = 'seismic.npy'
#woss(well_paths,wavelet_path,wavelet_size,seismic_path,ites=10000)
