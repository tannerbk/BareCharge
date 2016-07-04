#!/usr/bin/env python2.7

import numpy as np
import h5py
import ROOT
    
def read_hdf5(channel):
    samples = len(channel[0])
    waveforms = len(channel)
    dx = channel.attrs['horiz_interval']
    return samples, waveforms, dx

def draw_waveform(i, samples, dx, channel):

    canvas = ROOT.TCanvas()
    Waveform = ROOT.TH1F('waveform',';t (ns);ADC Counts',samples,0,samples*dx*1e9)
    for j in range(0, samples):
        adc = channel[i,j]
        Waveform.SetBinContent(j+1,adc)
        Waveform.SetBinError(j+1,1)
    print "Drawing waveform",
    print i
    yaxis = ("ACD Counts")
    Waveform.Draw("e")
    canvas.Update()
    Waveform.SetStats(0)
    Waveform.SetLineColor(ROOT.kViolet)

    raw_input()

def plot(filename, channel):

    f = h5py.File(filename)
    channel = f[channel]
    samples, waveforms, dx = read_hdf5(channel)

    for i in range(0, waveforms):
        draw_waveform(i, samples, dx, channel)

    return 0

if (__name__=="__main__"):
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    parser.add_argument('channel', help='Name of channel to process')
    args = parser.parse_args()
    plot(args.filename, args.channel)

