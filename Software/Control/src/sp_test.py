#!/usr/bin/env python
import sys, os
from ctypes import *
from sigproc import SigProc 
from ROOT import *
gROOT.LoadMacro("sp.C+")

from ROOT import SignalProcessor
from array import array
import re, time
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
# from process_signal import apply_wiener_filter
from scipy.signal import wiener
import cmath

def apply_wiener_filter(data,ich,n=16384):
    x = [data[ich*n+i] for i in range(n)]
    m = np.mean(x)
#     y = wiener([a-m for a in x], mysize=500)
    y = wiener([a-m for a in x], mysize=500, noise=0.0003)
    for i in range(n): data[ich*n+i] = y[i]+m

class bkgEstimator:
    def __init__(self):
        self.data = None
        self.npoints = None
        self.F = 59
        self.T = 16384/self.F
        self.get_data()
        self.scale = 0.6

    def get_data(self, evt=20):
        ch = TChain('tree1')
        ch.Add('data/fpgaLin/Feb09b_data_1881.root')
        n1 = ch.Draw('adc[19]','Entry$=={0:d}'.format(evt),'goff')
        v1 = ch.GetV1()
        x1 = np.array([v1[i] for i in range(n1)])

        n = int(n1/self.F)
        m1 = np.mean(x1)
        self.data = [np.mean(x1[i::n])-m1 for i in range(n)]
        self.npoints = n

        N = 16384
        self.Par = np.array([cmath.exp(-2*cmath.pi/N*59*k*1j) for k in range(N)])
        self.phase1 = int(self.get_phase(x1))

    def get_phase(self,x1):
        return np.angle(sum(x1*self.Par))/(2*cmath.pi)*self.T

    def correct(self, data,ich,n=16384): 
        phase = int(self.get_phase(data[ich*n:((ich+1)*n)]))
        print phase, self.phase1
#         for i in range(n): data[ich*n+i] -= self.data[(i-phase+self.phase1)%self.npoints]
        for i in range(n): data[ich*n+i] -= self.scale*self.data[(i+phase-self.phase1)%self.npoints]

    def show_data(self):
        plt.plot(self.data)

def test1():
    s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
    # data1 = s1.generate_adcDataBuf()
    # data2 = s1.generate_sdmDataBuf()
    # data1 = ((s1.ANALYSIS_WAVEFORM_BASE_TYPE * s1.nSamples) * s1.nAdcCh)()
    data1 = (s1.ANALYSIS_WAVEFORM_BASE_TYPE * (s1.nSamples * s1.nAdcCh))()
    # data1 = ((s1.ANALYSIS_WAVEFORM_BASE_TYPE * s1.nSamples) * s1.nAdcCh)()
    # print len(data1[0])
    print type(data1)
    # print type(byref(data1[0]))
    print type(pointer(data1))
    data1 = array('f',[0]*(16384*20))
    inRoot = 'data/fpgaLin/Feb09b_data_2.root'
    fout1 = TFile(inRoot,'read')
    tree1 = fout1.Get('tree1')
    tree1.SetBranchAddress('adc',data1)

    i = 37
    tree1.GetEntry(i)
    # print data1[3]

    sp1 = SignalProcessor()
    sp1.nAdcCh=20
    # sp1.sRanges.clear()
    # sp1.sRanges.push_back((0,3000))
    # sp1.sRanges.push_back((3000,8000))
    # sp1.sRanges.push_back((8000,15000))
    # sp1.measure_pulse(byref(data1))
    # sp1.measure_pulse(pointer(data1))
    # sp1.measure_pulse(data1)
    print sp1.nMeasParam

    f = 1000
    a = 16384*0.2*f*0.000001
    print 16384*0.2*f*0.000001
    n1 = int(1/(0.2*f*0.000001))
    print n1
    sp1.sRanges.clear()
    ip = 0
    while ip+n1<sp1.nSamples:
        iq = ip+n1
    #     if iq > sp1.nSamples: iq = sp1.nSamples
    #     a = (ip, iq)
        sp1.sRanges.push_back((ip, iq))
        ip = iq
    #
    sp1.measure_pulse(data1)
    for i in range(sp1.nAdcCh):
    #     print i, sp1.measParam[sp1.nMeasParam*i], sp1.measParam[sp1.nMeasParam*i+1], sp1.measParam[sp1.nMeasParam*i+2], sp1.measParam[sp1.nMeasParam*i+3]
        print i,
        for j in range(sp1.nMeasParam):
            print sp1.measParam[sp1.nMeasParam*i+j],
        print
    sp1.test2()


def test2():
    be1 = bkgEstimator()
#     be1.show_data()

    s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
    data1 = (s1.ANALYSIS_WAVEFORM_BASE_TYPE * (s1.nSamples * s1.nAdcCh))()
    data1 = array('f',[0]*(16384*20))


    inRoot = 'data/fpgaLin/Feb09b_data_1138.root'
    tag1 = 'Feb09b'
    if len(sys.argv)>1:
        if os.path.exists(sys.argv[1]):
            inRoot = sys.argv[1]
        elif os.path.exists('data/fpgaLin/Feb09b_data_'+sys.argv[1]+'.root'):
            inRoot = 'data/fpgaLin/'+tag1+'_data_'+sys.argv[1]+'.root'
        else:
            files = sorted([f for f in glob('data/fpgaLin/'+tag1+'_data_*.root')], key=lambda f:os.path.getmtime(f))

            a =  -1
            try:
                a = int(sys.argv[1])
            except TypeError:
                pass
            if time.time() - os.path.getmtime(files[-1]) < 10:
                print "dropping the latest file, which probably is still being written:", files[-1]
                if a!=0: files.pop()
                else: a = -1

            if abs(a)<len(files):
                inRoot = files[a]
            else:
                print "Index {0:d} out of range:{1:d}".format(a, len(files))
                return

    print "Using file:", inRoot
    fout1 = TFile(inRoot,'read')
    tree1 = fout1.Get('tree1')
    tree1.SetBranchAddress('adc',data1)

    print "Entries in the tree:", tree1.GetEntries()

    run = -1
    runPattern = '.*_data_(\d+).root'
    if runPattern is not None:
        m = re.match(runPattern, inRoot)
        if m:
            try:
                run = int(m.group(1))
            except ValueError:
                print "Run number not exatracted for file", iRoot

    i = 56
    ich = 19
    sp1 = SignalProcessor()
    sp1.fltParam.clear()
#     for x in [500, 500, 700, 2500]: sp1.fltParam.push_back(x)
#     for x in [500, 5, 15, 2500]: sp1.fltParam.push_back(x)
#     for x in [500, 50, 150, 2500]: sp1.fltParam.push_back(x)
#     for x in [30, 15, 50, 2500]: sp1.fltParam.push_back(x)
#     for x in [30, 10, 100, 2500]: sp1.fltParam.push_back(x)
#     for x in [30, 5, 100, 2500]: sp1.fltParam.push_back(x)
    for x in [30, 250, 350, 2500]: sp1.fltParam.push_back(x)
    sp1.x_thre = 0.05

    plt.ion()
    plt.show()
    fig, ax1 = plt.subplots(1, 1, figsize=(28, 10))
#     fig.set_size_inches(11,8)
    ax1.set_xlabel('time index')
    ax1.set_ylabel('U [V]', color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('U [V]', color='r')
    ax2.tick_params('y', colors='r')


#     for ievt in range(tree1.GetEntries()):

    NMax = tree1.GetEntries()
    ievt = 0
    while ievt< NMax:

        print "Event:", ievt
        tree1.GetEntry(ievt)

        va = data1[ich*sp1.nSamples:(ich+1)*sp1.nSamples]
        be1.correct(data1, ich)
#         apply_wiener_filter(data1, ich)

        sp1.measure_pulse2(data1, ich)

        vx = np.array([sp1.scrAry[i] for i in range(sp1.nSamples)])
        vo = data1[ich*sp1.nSamples:(ich+1)*sp1.nSamples]

        ax1.clear()
        ax2.clear()
        ax1.plot(va, label='Raw', color='b')
        ax1.plot(vo, label='Wiener', color='g')
        ax2.plot(vx, label='Filtered', color='r')
        ylim1 = ax1.get_ylim()
        ylim2 = ax2.get_ylim()

        x1 = min(ylim1[0], ylim2[0]+vo[0])
        x2 = max(ylim1[1], ylim2[1]+vo[0])
        print x1,x2
        ax1.set_ylim(x1,x2)
        ax2.set_ylim(x1-vo[0],x2-vo[0])

        print sp1.signals[ich].size()
        x1 = []
        y1 = []
        for s in sp1.signals[ich]:
            print s.idx, s.im, s.Q, s.w0, s.w1, s.w2
            x1.append(s.im)
            y1.append(s.Q)
            plt.axvline(x=s.im, linestyle='--', color='black')

        plt.text(0.04, 0.1, 'run {0:d} event {1:d}'.format(run, ievt), horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
        plt.xlim(auto=False)
        if x1: ax2.scatter(x1,y1, c="g", marker='o', s=220, label='Analysis')

        fig.tight_layout()
        plt.draw()
        plt.legend()
        plt.grid(True)
        plt.pause(0.001)


        while True:
            x = raw_input("Next:")
            if x=='q': sys.exit()
            elif len(x)>0 and x[0] == 's':
                for name in x.split()[1:]:
                    dirx = os.path.dirname(name)
                    if not os.path.exists(dirx): os.makedirs(dirx)
                    plt.savefig(name)
                    print "saved figure to", name
            else:
                try:
                    ievt = int(x)
                except ValueError:
                    ievt += 1
                break
            

#         fig.legend()
#         plt.show()

if __name__ == '__main__':
    test2()
