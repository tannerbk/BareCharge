#!/usr/bin/env python2.7

import h5py
import ROOT
    
def read_hdf5(channel):
    ''' 
    Read the hdf5 information
    '''

    samples = len(channel[0])
    dx = channel.attrs['horiz_interval']
    dy = channel.attrs['vertical_gain']
    return samples, dx, dy

def draw_waveform(i, samples, dx, dy, channel):
    '''
    Draw and individual scope waveform using ROOT
    '''

    canvas = ROOT.TCanvas()
    Waveform = ROOT.TH1F('waveform',';t (ns);Voltage(V)',samples,0,samples*dx*1e9)
    for j in range(0, samples):
        adc = channel[i,j]
        Waveform.SetBinContent(j+1,adc*dy)
        # One ADC-Count error bars
        Waveform.SetBinError(j+1,dy)
    Waveform.SetStats(0)
    Waveform.SetLineColor(ROOT.kBlack)
    Waveform.Draw("")
    canvas.Update()
    raw_input()

def plot(filename, channel):
    '''
    Open the hdf5 file and plot each waveform
    '''

    f = h5py.File(filename)
    channel = f[channel]
    samples, dx, dy = read_hdf5(channel)
    waveform_id = 0

    print "Drawing Waveforms, type 'CTRL-C' to exit."
    while True:
        print "Drawing waveform", waveform_id
        draw_waveform(waveform_id, samples, dx, dy, channel)
        waveform_id+=1

    return 0

if (__name__=="__main__"):
 
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    parser.add_argument('channel', help='Name of channel to process')
    args = parser.parse_args()
    plot(args.filename, args.channel)

