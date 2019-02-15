#!/usr/bin/env python
import sys, os, re
import matplotlib.pyplot as plt
from array import array
from multiprocessing import Pool
from glob import glob
from rootUtil import waitRootCmdX, useNxStyle
# from scipy.signal import find_peaks
from scipy.signal import find_peaks_cwt
import numpy as np
# from scipy import signal


from ROOT import *
gROOT.LoadMacro("sp.C+")

from ROOT import SignalProcessor


def check0():
    inRoot = sys.argv[1]      if len(sys.argv)>1 else "test.root"
    ich    = int(sys.argv[2]) if len(sys.argv)>2 else 0
    sp1 = SignalProcessor()
    sp1.nSamples = 16384 
    sp1.nAdcCh = 20
    freq = 1000
    n1 = int(1/(0.2*freq*0.000001))
    sp1.sRanges.clear()
    ip = 0
    dn = 2500%n1 ## 2500 is the expected position of the signal
    while ip+dn < sp1.nSamples:
        sp1.sRanges.push_back((ip, min(ip+n1, sp1.nSamples)))
        ip += n1

    data1 = array('f',[0]*(sp1.nSamples*sp1.nAdcCh))

    fout1 = TFile(inRoot,'read')
    tree1 = fout1.Get('tree1')
    tree1.SetBranchAddress('adc',data1)

    ### range
    eRange = None
    if len(sys.argv)>3:
        lx = []
        r = sys.argv[3].split(',')
        for ri in r:
            si = ri.find('-')
            if si == -1:
                lx.append(int(ri))
            else:
                end = tree1.GetEntries() if si==len(ri)-1 else int(ri[si+1:])
                lx += range(int(ri[:si]), end)
        eRange = lx
    else: eRange = range(tree1.GetEntries())

    ### plotting
    plt.ion()
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time (s)')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('exp', color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('sin', color='r')
    ax2.tick_params('y', colors='r')

#     for i in range(10):
    for i in eRange:
        print i
        tree1.GetEntry(i)
        ax1.clear()
        ax1.plot(data1[sp1.nSamples*ich:sp1.nSamples*(ich+1)])
        ax2.clear()

        sp1.fltParam.clear()
        for x in [500, 150, 200, 2500]: sp1.fltParam.push_back(x)
        sp1.measure_pulse(data1,ich)
        ax2.plot([sp1.scrAry[i] for i in range(sp1.nSamples)], 'r.')

        sp1.fltParam.clear()
        for x in [500, 15, 50, 2500]: sp1.fltParam.push_back(x)
#         for x in [500, 150, 200, 2500]: sp1.fltParam.push_back(x)
        sp1.measure_pulse(data1,ich)
#         ax2.plot([sp1.scrAry[i] for i in range(sp1.nSamples)], 'g.')
        
        vx = np.array([sp1.scrAry[i] for i in range(sp1.nSamples)])
        ax2.plot(vx, 'g.')

#         peakind = find_peaks_cwt(vx, np.arange(1,100))
#         peaks = find_peaks_cwt(vx, np.array([5, 10, 20,30,50,100,200,350]), min_snr=50)
#         print peaks, vx[peaks]
#         peaks, properties = signal.find_peaks(vx, prominence=1, width=20)
#         plt.plot(peaks, vx[peaks], "x")

        fig.tight_layout()
        plt.draw()
        plt.grid(True)
        plt.pause(0.001)
        x = raw_input("Press [enter] to continue.")
        if x=='q': break


def check1(inRoot, chRange=None, nEvt=None, oTag='_tt'):
    '''Use various filters to extract the results'''

    if os.path.exists(inRoot.rstrip('.root')+oTag+'.root'):
        print "has:", inRoot
        return

#     print inRoot
#     return
    sp1 = SignalProcessor()
    sp1.nSamples = 16384 
    sp1.nAdcCh = 20
    freq = 1000
    n1 = int(1/(0.2*freq*0.000001))
    sp1.sRanges.clear()
    ip = 0
    dn = 2500%n1 ## 2500 is the expected position of the signal
    while ip+dn < sp1.nSamples:
        sp1.sRanges.push_back((ip, min(ip+n1, sp1.nSamples)))
        ip += n1

    data1 = array('f',[0]*(sp1.nSamples*sp1.nAdcCh))

    fin1 = TFile(inRoot,'read')
    tree1 = fin1.Get('tree1')
    tree1.SetBranchAddress('adc',data1)
    if nEvt is None: nEvt = tree1.GetEntries()
    if chRange is None: chRange = range(sp1.nAdcCh)

    fout1 = TFile(inRoot.rstrip('.root')+oTag+'.root','recreate')
    tup1 = TNtuple('tup1',"filter analysis tuple",'evt:fR:fW:ich:b:bE:im:idx:A')

    ### start processing
    INTV = 1000 if nEvt>10000 else max(nEvt/10, 1)
    for ievt in range(nEvt):
        tree1.GetEntry(ievt)

        if ievt%INTV==0: print ievt, ' events processed'

        for R in range(50,300, 50):
            for W in range(50, 600, 50):
                sp1.fltParam.clear()
                for x in [500, R, W, -1.]: sp1.fltParam.push_back(x)
                sp1.measure_pulse(data1)

                for ich in chRange:
                    itmp = sp1.nMeasParam*ich
                    for j in range(sp1.nMeasParam/2-1):
                        tup1.Fill(ievt, R, W, ich, sp1.measParam[itmp], sp1.measParam[itmp+1], j, sp1.measParam[itmp+2*j+2], sp1.measParam[itmp+2*j+3])
    tup1.Write()

    ### save the ENC results
    tup2 = TNtuple('tup2',"filter analysis tuple enc",'ich:fR:fW:m:mE:sigma:sigmaE:fq:fStatus')
    for R in range(50,300, 50):
        for W in range(50, 600, 50):
            for ich in chRange:
                gDirectory.Delete('h1*')
                tup1.Draw('A>>h1','int(fW)=={0:d}&&int(fR)=={1:d}&&ich=={2:d}'.format(W,R,ich),'goff')
                h1 = gDirectory.Get('h1')
                r = h1.Fit('gaus', 'S0')
                fun1 = h1.GetFunction('gaus')
                tup2.Fill(ich, R, W, fun1.GetParameter(1),fun1.GetParError(1), fun1.GetParameter(2),fun1.GetParError(2), r.Prob(), r.Status())

    tup2.Write()
    fout1.Close()

def run1():
    inRoot = sys.argv[1]      if len(sys.argv)>1 else "test.root"
    ich    = sys.argv[2].split(',') if len(sys.argv)>2 else None
    nEvt   = int(sys.argv[3]) if len(sys.argv)>3 else None

    check1(inRoot, ich, nEvt)

def process_check1():
    p = Pool(6)
    p.map(check1, glob('data/fpgaLin/Jan25a_C2_*_f1000.root'))

def checkENC():
    ich = 0
    fin1 = TFile('out_filter_check.root','read')
    tup1 = fin1.Get('tup1')

    fout1 = TFile('out_filter_check_enc.root','recreate')
    tup2 = TNtuple('tup2',"filter analysis tuple enc",'ich:fR:fW:m:mE:sigma:sigmaE:fq:fStatus')

    for R in range(50,300, 50):
        for W in range(50, 600, 50):
            gDirectory.Delete('h1*')
            tup1.Draw('A>>h1','int(fW)=={0:d}&&int(fR)=={1:d}'.format(W,R),'goff')
            h1 = gDirectory.Get('h1')
#             print R, W, h1.GetEntries()
            r = h1.Fit('gaus', 'S0')
            fun1 = h1.GetFunction('gaus')
            tup2.Fill(ich, R, W, fun1.GetParameter(1),fun1.GetParError(1), fun1.GetParameter(2),fun1.GetParError(2), r.Prob(), r.Status())
#             h1.Draw()
#             gPad.Update()
#             a = raw_input()
#             if a == 'q': return
#             waitRootCmdX()

    tup2.Write()
    fout1.Close()

def checkENC_d():
    ich = 0
    fin1 = TFile('data/fpgaLin/Jan25a_C2_200mV_f1000_tt.root','read')
    tup1 = fin1.Get('tup1')

    for R in range(50,300, 50):
        for W in range(50, 600, 50):
            gDirectory.Delete('h1*')
            tup1.Draw('A>>h1','int(fW)=={0:d}&&int(fR)=={1:d}&&ich=={2:d}'.format(W,R,ich),'goff')
            h1 = gDirectory.Get('h1')
#             print R, W, h1.GetEntries()
            r = h1.Fit('gaus', 'S0')
            fun1 = h1.GetFunction('gaus')
#             tup2.Fill(ich, R, W, fun1.GetParameter(1),fun1.GetParError(1), fun1.GetParameter(2),fun1.GetParError(2), r.Prob(), r.Status())
            if r.Prob()<0.05:
                h1.Draw()
                gPad.Update()
                a = raw_input()
                if a == 'q': return
#             waitRootCmdX()

#     tup2.Write()
#     fout1.Close()
def scan_run():
    with open('flt_test.ttl','w') as fout1:
        fout1.write('ich/I:Vin/F:W/I:R/I:enc')

        for fname in glob('data/fpgaLin/Jan25a_C2_*mV_f1000_tt.root'):
            fin1 = TFile(fname,'read')
            tup2 = fin1.Get('tup2')
            if not tup2:
                print fname, 'tree'
                continue

            m = re.match(r".*_C2_(\d+)mV_f1000_tt.root", fname)
            if m is None:
                print fname
                continue
            vin = m.group(1)

            for ich in range(20):
                n = tup2.Draw("sigma/m:fW:fR","ich=={0:d}".format(ich),"goff")
                if n == 0:
                    print fname, n, ich
                    continue
                enc = tup2.GetV1()
                fW = tup2.GetV2()
                fR = tup2.GetV3()

                a = TMath.LocMin(n, tup2.GetV1())
                fout1.write('\n'+' '.join([str(ich),vin,str(fW[a]), str(fR[a]), str(enc[a])]))

def scan_test():
    ich = 0
    fin1 = TFile('data/fpgaLin/Jan25a_C2_100mV_f1000_tt.root','read')
    tup2 = fin1.Get('tup2')

    n = tup2.Draw("sigma/m:fW:fR","ich==2","goff")
    enc = tup2.GetV1()
    fW = tup2.GetV2()
    fR = tup2.GetV3()

    a = TMath.LocMin(n, tup2.GetV1())
    print a, enc[a], fW[a], fR[a]
    for i in range(n): print '--', i, enc[i]

def main():
    inRoot = sys.argv[1]      if len(sys.argv)>1 else "test.root"
    ich    = int(sys.argv[2]) if len(sys.argv)>2 else 0
    sp1 = SignalProcessor()
    sp1.nSamples = 16384 
    sp1.nAdcCh = 20
    sp1.fltParam.clear()
    for x in [50, 150, 200, -1.]: sp1.fltParam.push_back(x)

    freq = 1000
    n1 = int(1/(0.2*freq*0.000001))
    sp1.sRanges.clear()
    ip = 0
    dn = 2500%n1 ## 2500 is the expected position of the signal
    while ip+dn < sp1.nSamples:
        sp1.sRanges.push_back((ip, min(ip+n1, sp1.nSamples)))
        ip += n1

    data1 = array('f',[0]*(sp1.nSamples*sp1.nAdcCh))

    fout1 = TFile(inRoot,'read')
    tree1 = fout1.Get('tree1')
    tree1.SetBranchAddress('adc',data1)

    ### range
    eRange = None
    if len(sys.argv)>3:
        lx = []
        r = sys.argv[3].split(',')
        for ri in r:
            si = ri.find('-')
            if si == -1:
                lx.append(int(ri))
            else:
                end = tree1.GetEntries() if si==len(ri)-1 else int(ri[si+1:])
                lx += range(int(ri[:si]), end)
        eRange = lx
    else: eRange = range(tree1.GetEntries())

    ### plotting
    plt.ion()
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('time (s)')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('exp', color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('sin', color='r')
    ax2.tick_params('y', colors='r')

#     for i in range(10):
    for i in eRange:
        print i
        tree1.GetEntry(i)
        sp1.measure_pulse(data1,ich)

        ax1.clear()
        ax1.plot(data1[sp1.nSamples*ich:sp1.nSamples*(ich+1)])
        ax2.clear()
        ax2.plot([sp1.scrAry[i] for i in range(sp1.nSamples)], 'r.')
        fig.tight_layout()
        plt.draw()
        plt.pause(0.001)
        x = raw_input("Press [enter] to continue.")
        if x=='q': break

if __name__ == '__main__':
#     main()
#     scan_test()
#     scan_run()
    check0()
#     run1()
#     process_check1()
#     checkENC()
#     checkENC_d()
