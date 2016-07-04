#!/usr/bin/env python2.7

import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import colorsys
from matplotlib.patches import Wedge 
from matplotlib.patches import Rectangle

termination_ohms = 50

def read_hdf5(channel):
    samples = len(channel[0])
    dx = channel.attrs['horiz_interval']
    dy = channel.attrs['vertical_gain']
    return samples, dx, dy

def draw_waveform(i, samples, dx, dy, channel, chn):

    pedestal = np.mean(channel[i,:1000])
    vmin = np.min(channel[i,:])
    q = 0
    for j in range(1000, samples):
        v = (channel[i,j] - pedestal)*dy
        q += (v*((-1000*dx*1e9)/termination_ohms))
  
    print "Waveform:",
    print i
    charge_geometry(i, q, chn)
    #time_geometry(i,vmin,chn)

def plot(filename):

    f = h5py.File(filename)
    for k in range(0, 50):
        print ' '
        chn = 0
        for name in ['channel%i' % i for i in range(1,5)]:
            chn+=1
            if name in f:
                print name,
                print '     ',
                channel = f[name]
                samples, dx, dy = read_hdf5(channel)
                draw_waveform(k, samples, dx, dy, channel, chn)

    return 0

def charge_geometry(waveform_number, q, channel_number):
    
    Color = q/(50*1.6)
 
    # probably should use Colormap.Normalize
    if q/(50*1.6) > 1.0:
        Color = 1.0
    if q < 0.0:
        Color = 0.0

    r, g, b, a = mpl.cm.jet(Color)

    # Draw at start
    if waveform_number == 0 and channel_number == 1:
        fig = plt.figure(figsize=(10,10))
        ax1 = fig.add_axes([0.05, 0.95, 0.9, 0.15])
        cmap = mpl.cm.jet
        norm = colors.Normalize(vmin = 0, vmax = 35)
        cb1 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,norm=norm,orientation='horizontal')
        cb1.set_label('Number of Photons')
        ax2 = fig.add_axes([0.05, 0.1, 0.9, 0.75])
        scint=plt.Circle((.5,.5),.2,color='b')
        plt.axis('off')
        fig.gca().add_artist(scint)
    else:  
        fig = plt.gcf()
    
    print "PE:",
    print q/1.6
    
    # channels displayed as:
    #       ch2
    #    ch3 O ch1
    #       ch4 
    if channel_number == 1:
        w1 = Wedge((0.82,0.5),0.12,90,270,color=(r,g,b))
        r1 = Rectangle((0.82,0.45),0.1,0.1,color=(r,g,b))
        #pmt1=plt.Circle((0.82,0.5),.12,color=(r,g,b))
        fig.gca().add_artist(w1) 
        fig.gca().add_artist(r1) 
    if channel_number == 2:
        w2 = Wedge((0.5,0.82),0.12,180,0,color=(r,g,b))
        r2 = Rectangle((0.45,0.82),0.1,0.1,color=(r,g,b))
        #pmt2=plt.Circle((0.5,0.82),.12,color=(r,g,b))
        fig.gca().add_artist(w2)
        fig.gca().add_artist(r2)
    if channel_number == 3:
        w3 = Wedge((0.18,0.5),0.12,270,90,color=(r,g,b)) 
        r3 = Rectangle((0.08,0.45),0.1,0.1,color=(r,g,b))
        #pmt3=plt.Circle((0.18,0.5),.12,color=(r,g,b))
        fig.gca().add_artist(w3)
        fig.gca().add_artist(r3)
    if channel_number == 4:
        w4 = Wedge((0.5,0.18),0.12,0,180,color=(r,g,b))
        r4 = Rectangle((0.45,0.08),0.1,0.1,color=(r,g,b))
        #pmt4=plt.Circle((0.5,0.18),.12,color=(r,g,b))
        fig.gca().add_artist(w4)
        fig.gca().add_artist(r4)
        fig.show()
        print 'Enter to continue.',
        print 'CTRL-C+Enter to exit.'
        raw_input()

if (__name__=="__main__"):
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    args = parser.parse_args()
    plot(args.filename)

