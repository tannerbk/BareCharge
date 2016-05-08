#!/usr/bin/env python2.7

import numpy as np
import h5py
import ROOT
    

'''  Get the waveform attributes
   : param filename:   data file
   : param channel:    data channel
     Returns number of samples, number of waveforms
     and horizontal resolution
'''
def read_hdf5(channel):
    samples = len(channel[0])
    waveforms = len(channel)
    dx = channel.attrs['horiz_interval']
    return samples, waveforms, dx

'''  Draw waveforms with the debug option 
   : param i:            waveform number 
   : param samples:      number of samples
   : param channel:      channel
     Returns TH1 Waveform of scope trace
'''

def draw_waveform(i, samples, dx, channel):

    Waveform = ROOT.TH1F('waveform',';t (ns);ADC Counts',samples,0,samples*dx*1e9)
    for j in range(0, samples):
        adc = channel[i,j]
        Waveform.SetBinContent(j+1,adc)
        Waveform.SetBinError(j+1,1)
    print "Drawing waveform",
    print i
    yaxis = ("ACD Counts")
    Waveform.Draw("e")
    Waveform.SetStats(0)
    Waveform.SetLineColor(ROOT.kViolet)
    canvas.Update()

    raw_input()

''' In debug mode this calls draw_waveform, which will continue to
    print waveforms to the screen. Cancel with CTRL-C.
'''
def plot(filename, channel):

    f = h5py.File(filename)
    channel = f[channel]
    samples, waveforms, dx = read_hdf5(channel)

    global canvas
    canvas = ROOT.TCanvas()

    for i in range(0, waveforms):
        draw_waveform(i, samples, dx, channel)

    return 0

''' Input arguments to process
   : arg filename:        data file to process
   : arg channel:         channel to process
'''
if (__name__=="__main__"):
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    parser.add_argument('channel', help='Name of channel to process')
    args = parser.parse_args()
    plot(args.filename, args.channel)

