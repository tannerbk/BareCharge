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
        y = np.exp(-0.5 * (np.log(x[0] / par[0]) / par[1])**2.0)
        y *= par[2] / (x[0] * par[1] * np.sqrt(2.0 * np.pi))
        p = np.exp(-0.5 * (np.log(x[0] / par[3]) / par[4])**2.0)
        p *= par[5] / (x[0] * par[4] * np.sqrt(2.0 * np.pi))
        q = y + p
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
   : param i:            waveform number 
   : param samples:      number of samples
   : param channel:      channel
   : param dy:           veritcal resolution
   : param min_bin:      bin correspong to peak
   : param pedestal:     length of pedestal window in time bins
   : param num_channels: number of channels to draw (max 2)
   : param fit:          fit the trigger waveform w/ lognormal
     Returns TH1 Waveform of scope trace
'''

def draw_waveform(f, i, samples, channel, dy, dx, min_bin, pedestal, num_channels, fit):
    global canvas
    canvas = ROOT.TCanvas()
    pad1 = ROOT.TPad("pad1","pad1",0,0.5,1,1)
    pad2 = ROOT.TPad("pad2","pad2",0,0.5,1,0.0)
    Waveform = ROOT.TH1F('waveform',';t(ns);Voltage',samples,0,samples*0.1)
    channel_list = channels(num_channels)
    for ch in channel_list:
        c = f[ch]
        for j in range(0, samples):
            voltage = c[i,j]*dy - pedestal
            Waveform.SetBinContent(j+1,voltage)
            Waveform.SetBinError(j+1,dy)
        print "Drawing waveform",
        print i
        if ch == "channel1":
            print ch
            yaxis = ("Voltage (V) "+ch)
            pad1.SetBottomMargin(0)
            pad1.Draw()
            pad1.cd()
            Waveform.Draw("e")
            Waveform.SetStats(0)
            Waveform.SetLineColor(1)
            Waveform.GetYaxis().SetTitleOffset(1.20)
            Waveform.GetYaxis().SetTitle(yaxis)
            # do the fit for channel 1, and draw
            if fit==True:
                tmin, bottom_range, top_range = min_bin*dx*1e9, min_bin*dx*1e9-3.0, min_bin*dx*1e9+4.0
                fit = ROOT.TF1("lognormal",lognormal(),bottom_range,top_range,6)
                fit.SetParameter(0, tmin)
                fit.SetParameter(1, 0.02)
                fit.SetParameter(2, 0.1)
                fit.SetParameter(3, tmin)
                fit.SetParameter(4, 0.02)
                fit.SetParameter(5, 0.1)
                fit.SetParLimits(0, tmin-10, tmin+10)
                Waveform.Fit(fit,"0")
                fit.SetLineColor(ROOT.kSpring)
                fit.Draw("same")
                Waveform.GetXaxis().SetRangeUser(bottom_range-10,top_range+10)
            canvas.Update()
        else:
            print ch
            yaxis = ("Voltage (V) "+ch)
            canvas.cd()
            pad2.SetTopMargin(0)
            pad2.SetBottomMargin(0.1)
            pad2.Draw()
            pad2.cd()
            Waveform.Draw("e same")
            Waveform.SetStats(0)
            Waveform.SetLineColor(2)
            Waveform.GetYaxis().SetTitleOffset(1.20)
            Waveform.GetYaxis().SetTitle(yaxis)
            Waveform.GetXaxis().SetRangeUser(0,500)
            canvas.Update()
    #canvas.Update()
    raw_input()
    pad1.Close()
    pad2.Close()

''' In debug mode this calls draw_waveform, which will continue to
    print waveforms to the screen. Cancel with CTRL-C.
    In normal mode, this builds charge and pedestal histograms
    and write thems to a file.
'''
def plot(filename, channel, pedestal_window, num_channels, debug=True, fit=False):
    Tf = ROOT.TFile("output.root", 'recreate')
    Charge = ROOT.TH1F('charge',';Charge(pC);Counts', 1250,-1,10)
    Pedestal = ROOT.TH1F('pedestal',';Voltage(V);Counts', 50000,-1,1)

    f = h5py.File(filename)
    channel = f[channel]
    print "Number of channels", 
    print num_channels
    samples, waveforms, dy, dx = read_hdf5(f, channel)
        
    for i in range(0,waveforms):
        charge = 0
        min_bin, pedestal = peak_finder(i, samples, channel, dy, pedestal_window)
        if debug:
            draw_waveform(f, i, samples, channel, dy, dx, min_bin, pedestal, num_channels, fit)
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
   : arg --num_channels:  number of channels to draw (currently only works for 2)
   : arg --fit:           fit the trigger signal with a lognormal
'''
if (__name__=="__main__"):
    import argparse

    parser = argparse.ArgumentParser(description='Draw waveforms')
    parser.add_argument('filename', help='Name of HDF5 file to process')
    parser.add_argument('channel', help='Name of channel to process')
    parser.add_argument('pedestal_window', type=int, help='Length of pedestal_window')
    parser.add_argument('--debug', action='store_true', help='Show debugging plots')
    parser.add_argument('--num_channels', type=int, default=1,
                         help='Number of channels to process')
    parser.add_argument('--fit', type=bool, default=False)
    args = parser.parse_args()
    plot(args.filename, args.channel, args.pedestal_window, args.num_channels, args.debug, args.fit)

