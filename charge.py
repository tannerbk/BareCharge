#!/usr/bin/env python2.7
import numpy as np
import h5py
import ROOT
import time

''' Code to analyze .hdf5 files 
    that store the fully digitized
    scope waveform 
'''
# TO DO: improve speed

canvas = None
termination_ohms = 50 

# Lognormal distribution for the pulse fit 
class lognormal:
    def __call__(self,x,par):
        q = np.exp(-0.5 * (np.log(x[0] / par[0]) / par[1])**2.0)
        q *= par[2] / (x[0] * par[1] * np.sqrt(2.0 * np.pi))
        return -q

'''  Get the waveform attributes
   : param filename:   data file
   : param channel:    data channel
     Returns number of samples, number of waveforms
     horizontal, and verical resolutions
'''
def read_hdf5(filename, channel):
    samples = len(channel[0])
    waveforms = len(channel)
    dy = channel.attrs['vertical_gain']
    dx = channel.attrs['horiz_interval']
    return samples, waveforms, dy, dx

'''  Find the bin corresponding to the minimum voltage 
     Also returns the pedestal for the given waveform 
'''
def peak_finder(i, samples, channel, dy, pedestal_window):
    min_voltage = 0.0
    pedestal = np.mean(channel[i,:pedestal_window])*dy
    for j in range(0, samples):
        voltage = channel[i,j]*dy - pedestal
        if voltage < min_voltage:
            min_voltage = voltage
            min_bin = j

    return min_bin, pedestal

''' Returns an array with all the channels
    Ex: ['channel1','channel2',channel3']
'''
def channels(num_channels):
    channel = []
    for i in range(1,num_channels+1):
        channel.append("channel"+str(i))

    return channel

'''  Draw waveforms with the debug option 
   : param i:          waveform number 
   : param samples:    number of samples
   : param channel:    channel
   : param dy:         veritcal resolution
   : param min_bin:    bin correspong to peak
     Returns TH1 Waveform of scope trace
'''
def draw_waveform(i, samples, channel, dy, dx, min_bin, pedestal):
    global canvas
    canvas = ROOT.TCanvas();
    Waveform = ROOT.TH1F('waveform',';t(ns);Voltage',samples,0,samples*0.1)
    for j in range(0, samples):
        voltage = channel[i,j]*dy - pedestal
        Waveform.SetBinContent(j+1,voltage)
        Waveform.SetBinError(j+1,dy)
    print "Drawing waveform",
    print i
    Waveform.Draw("e")
    Waveform.SetLineColor(1)
    tmin, bottom_range, top_range = min_bin*dx*1e9, min_bin*dx*1e9-15.0, min_bin*dx*1e9+5.0

    # do the fit, and draw
    f = ROOT.TF1("lognormal",lognormal(),bottom_range,top_range,3)
    f.SetParameter(0, tmin)
    f.SetParameter(1, 0.02)
    f.SetParameter(2, 0.1)
    f.SetParLimits(0, tmin-10, tmin+10)
    Waveform.Fit(f,"0")
    f.SetLineColor(ROOT.kSpring)
    f.Draw("same")

    canvas.Update()
    raw_input()

''' In debug mode this calls draw_waveform, which will continue to
    print waveforms to the screen. Cancel with CTRL-C.
    In normal mode, this builds charge and pedestal histograms
    and write thems to a file.
'''
def plot(filename, channel, pedestal_window, num_channels, debug=True):
    Tf = ROOT.TFile("output.root", 'recreate')
    Charge = ROOT.TH1F('charge',';Charge(pC);Counts', 1250,-1,10)
    Pedestal = ROOT.TH1F('pedestal',';Voltage(V);Counts', 50000,-1,1)

    f = h5py.File(filename)
    channel = f[channel]
    print "Number of channels", 
    print num_channels
    samples, waveforms, dy, dx = read_hdf5(f, channel)
        

    print "Analyzing waveform "
    for i in range(0,waveforms):
        print i,
        charge = 0
        min_bin, pedestal = peak_finder(i, samples, channel, dy, pedestal_window)
        if debug:
            draw_waveform(i, samples, channel, dy, dx, min_bin, pedestal)
        for j in range(min_bin-250,min_bin+250):
            voltage = channel[i,j]*dy - pedestal
            charge += (voltage*((-1000*dx*1e9)/termination_ohms))
        Charge.Fill(charge)
        Pedestal.Fill(pedestal)

    Tf.cd()
    Charge.Write()
    Pedestal.Write()
    Tf.Close()
    return 0

''' Input arguments to process
   : arg filename:        data file to process
   : arg channel:         channel to process
   : arg pedestal_window: length of the pedestal window, in time-bins
   : arg --debug:         write waveforms to screen
'''
if (__name__=="__main__"):
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    parser.add_argument('channel', help='Name of channel to process')
    parser.add_argument('pedestal_window', type=int, help='Length of pedestal_window')
    parser.add_argument('--debug', action='store_true', help='Show debugging plots')
    parser.add_argument('--num_channels', type=int, help='Number of channels to process')
    args = parser.parse_args()
    num_channels = 1
    plot(args.filename, args.channel, args.pedestal_window, num_channels, args.debug)

