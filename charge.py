#!/usr/bin/env python2.7
import numpy as np
import h5py
import ROOT

canvas = None

def read_hdf5(filename, channel,  debug=True):
    f = h5py.File(filename)
    channel = f[channel]
    
    samples = len(channel[0]) # samples per waveform 
    waveforms = len(channel) # number of waveforms 
    termination_ohms = 50 
    dy = channel.attrs['vertical_gain'] # vertical resolution 
    dx = channel.attrs['horiz_interval'] # time-bin

    global canvas
    if debug and canvas is None:
        canvas = ROOT.TCanvas();

    Tf = ROOT.TFile("output.root", 'recreate')
    Waveform = ROOT.TH1F('waveform',';t(ns);Voltage',samples,0,samples*0.1)
    Charge = ROOT.TH1F('charge',';Charge(pC);Counts', 1250,-1,10)

    voltage = 0
    for i in range(0,waveforms):
        charge = 0
        pedestal = np.mean(channel[i,:])*dy
        # Builds waveform 
        for j in range(0,samples):
            voltage = channel[i,j]*dy - pedestal
            Waveform.SetBinContent(j+1,voltage)
            Waveform.SetBinError(j+1,dy)
        # Hard coded signal window, use to intergrate PMT signal 
        for j in range(2800,3300):
            voltage = channel[i,j]*dy - pedestal
            charge = charge+(voltage*((-1000*dx*1e9)/termination_ohms))
        Charge.Fill(charge)
        # Debug mode prints waveforms to screen 
        if debug:
             print "Waveform",
             print i 
             Waveform.Draw("e")
             Waveform.SetLineColor(1)
             canvas.Update()
             raw_input()

    Tf.cd()
    Waveform.Write()
    Charge.Write()
    Tf.Close()
    return 0

if (__name__=="__main__"):
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms.')
    parser.add_argument('--debug', action='store_true',
                        help='Show debugging plots')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    parser.add_argument('channel', help='Name of channel to process')
    args = parser.parse_args()
    read_hdf5(args.filename, args.channel, args.debug)
