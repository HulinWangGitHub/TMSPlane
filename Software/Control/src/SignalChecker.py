#!/usr/bin/env python
import socket, os
# import argparse
# import pprint
from command import *
from sigproc import SigProc 
import time
import array
import glob
from ROOT import *
from subprocess import call
from math import isnan

def waitRootCmdY():
    a = raw_input("waiting...")

def useLHCbStyle0():
    pass

try:
    from rootUtil import waitRootCmdX, useLHCbStyle
except ImportError:
    waitRooCmdX = waitRootCmdY
    useLHCbStyle = useLHCbStyle0 


class SignalChecker:
    def __init__(self):
        self.control_ip_port = "192.168.2.3:1024"
        self.cmd = Cmd()
        self.s = None
        self.connected = False
        self.fileSuffix = '.1'
#         self.connect()

    def connect(self):
        ctrlipport = self.control_ip_port.split(':')
        self.s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        self.s.connect((ctrlipport[0],int(ctrlipport[1])))
        self.connected = True

    def take_samples(self, n=10, name="sample_{0:d}"):
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        # FPGA internal fifo : 512 bit width x 16384 depth
        nWords = (512 // 32) * 16384
        nBytes = nWords * 4
        buf = bytearray(nBytes)

        for i in range(n):
            self.s.sendall(self.cmd.send_pulse(1<<2));
            time.sleep(0.5)

            name1 = name.format(i)
            ret = self.cmd.acquire_from_datafifo(self.s, nWords, buf)
            s1.demux_fifodata(ret,data1,data2)
            s1.save_data([name1+'.adc', name1+'.sdm'], data1, data2)

    def take_samples2(self, n=10, outRootName='test_sample.root', runNumber=0):
        if not self.connected: self.connect()
        nMonitor = 20

        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        T = array.array('i',[0])
        V = array.array('i',[0])
        if self.fileSuffix:
            while os.path.exists(outRootName): outRootName += self.fileSuffix
        fout1 = TFile(outRootName,'recreate')
        tree1 = TTree('tree1',"data: {0:d} channel, {1:d} samples".format(s1.nAdcCh, s1.nSamples))
        tree1.Branch('T',T,'T/i')
        tree1.Branch('V',V,'V/I')
        tree1.Branch('adc',data1, "adc[{0:d}][{1:d}]/F".format(s1.nAdcCh, s1.nSamples))

        # FPGA internal fifo : 512 bit width x 16384 depth
        nWords = (512 // 32) * 16384
        nBytes = nWords * 4
        buf = bytearray(nBytes)

        status = 0
        try:
            for i in range(n):
                if i%100==0:
                    print(str(i)+' samples taken.')
                    try:
                        with open('/home/dlzhang/work/repos/TMSPlane2/Software/Control/src/.pulse_status') as f1:
                            V[0] = int(f1.readline().rstrip())
                    except:
                        V[0] = -999
                self.s.sendall(self.cmd.send_pulse(1<<2));

                T[0] = int(time.time())
                ret = self.cmd.acquire_from_datafifo(self.s, nWords, buf)
                s1.demux_fifodata(ret,data1,data2)
                tree1.Fill()

                if i%nMonitor == 1: tree1.AutoSave("SaveSelf");
        except KeyboardInterrupt:
            status = 1

        fout1.Write()
        return status

    def show_signal(self):
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        # FPGA internal fifo : 512 bit width x 16384 depth
        nWords = (512 // 32) * 16384
        nBytes = nWords * 4
        buf = bytearray(nBytes)


        self.s.sendall(self.cmd.send_pulse(1<<2));
        time.sleep(0.5)

        ret = self.cmd.acquire_from_datafifo(self.s, nWords, buf)
        s1.demux_fifodata(ret,data1,data2)

        mx = s1.measure_pulse(data1, fltParam=[500, 150, 200, -1.])
        for x in mx:
            for i in range(len(x)):
                print x[i]
        d = [x[1] for x in mx]
        self.plot_data(d)

    def plot_data(self,d):
        if len(d)<19:
            print('Need 19 entries')
            return
        from ROOT import TH2Poly, gStyle
        gStyle.SetOptStat(0)

        cv1 = TCanvas('cv1','cv1',600,500)
        cv1.SetRightMargin(0.2)
        hc = TH2Poly();
        hc.SetTitle("TMS19Plane");
        hc.Honeycomb(-4.3,-4,1,5,5);
        listA = [(0,0),(2,0),(1,1.5),(-1,1.5),(-2,0),(-1,-1.5),(1,-1.5),(4,0),(3,2),(1,3),(0,3),(-1,3),(-3,2),(-4,0),(-3,-2),(-1,-3),(0,-3),(1,-3),(3,-2)]

        for i in range(len(listA)):
            hc.Fill(listA[i][0],listA[i][1],d[i])

        hc.Draw("text colz0");
        waitRootCmdX()

    def check_file(self, fname):
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        s1.read_data([fname,''], data1, data2)

        mx = s1.measure_pulse(data1, fltParam=[500, 150, 200, -1.])
        d = [x[2] for x in mx]
        self.plot_data(d)


    def check_event(self, rootfile, idxs):
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()

        fout1 = TFile(rootfile,'read')
        tree1 = fout1.Get('tree1')
        tree1.SetBranchAddress('adc',data1)
        for i in idxs:
            tree1.GetEntry(i)
            mx = s1.measure_pulse(data1, fltParam=[500, 150, 200, -1.])
            d = [x[3] for x in mx]
            self.plot_data(d)

            waitRootCmdX()



    def show_sample(self, fnameP='/data/Samples/TMSPlane/Dec26/sample_{0:d}.adc', Ns=10, ich=8):
        from ROOT import TGraph, gPad, TLatex
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        data3 = ((s1.ANALYSIS_WAVEFORM_BASE_TYPE * s1.nSamples)*Ns)()
        dataList = []
        for i in range(Ns):
            fname = fnameP.format(i)
            s1.read_data([fname,''], data1, data2)
            x = array.array('f', range(s1.nSamples))
#             s1.filters_trapezoidal(data1[ich], data3[i], [100,100,200,-1])
            s1.filters_trapezoidal(data1[ich], data3[i], [100,100,100,0.000722])
            g1 = TGraph(s1.nSamples, x, data3[i])
#             g1 = TGraph(s1.nSamples, x, data1[ich])

            opt = 'AP' if i==0 else "PSame"
            g1.Draw(opt+' PLC PMC')
            dataList.append(g1)

        lt = TLatex()
        lt.DrawLatexNDC(0.2,0.8,"Chan={0:d}".format(ich))
        gPad.Update()
        waitRootCmdX()

    def check_enc(self, filePattern):
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        files = glob.glob(filePattern)
        with open('test1bbb.dat','w') as fout:
            fout.write(':'.join(['sID/I', 'ch/I', 'B/F','dB/F','idx/I','A/F']))
            for fname in files:
                s1.read_data([fname,''], data1, data2)
                mxx = s1.measure_pulse(data1, fltParam=[500, 150, 200, -1.])
                for ch in range(s1.nAdcCh):
                    mx = mxx[ch]
                    fout.write('\n'+' '.join([fname[fname.find('_')+1:-4],str(ch), str(mx[0]), str(mx[1]), '{0:d}'.format(int(mx[2])), str(mx[3])]))

    def check_enc2(self, inRoot, outText):
        s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
        data1 = s1.generate_adcDataBuf()
        data2 = s1.generate_sdmDataBuf()

        fout1 = TFile(inRoot,'read')
        tree1 = fout1.Get('tree1')
#         tree1.SetBranchAddress('ch0',data1)
        tree1.SetBranchAddress('adc',data1)
        with open(outText,'w') as fout:
            fout.write(':'.join(['sID/I', 'ch/I', 'B/F','dB/F','idx/I','A/F']))

            for i in range(tree1.GetEntries()):
                if i%100==0: print(str(i)+' entries processed')
                tree1.GetEntry(i)
                mxx = s1.measure_pulse(data1, fltParam=[500, 150, 200, -1.])
                for ch in range(s1.nAdcCh):
                    mx = mxx[ch]
                    if isnan(mx[1]): print("xxxx")
                    fout.write('\n'+' '.join([str(i),str(ch), str(mx[0]), str(mx[1]), '{0:d}'.format(int(mx[2])), str(mx[3])]))

def text2root(spattern, irange, outname):
    s1 = SigProc(nSamples=16384, nAdcCh=20, nSdmCh=19, adcSdmCycRatio=5)
    data1 = s1.generate_adcDataBuf()
    data2 = s1.generate_sdmDataBuf()

    fout1 = TFile(outname,'recreate')
    tree1 = TTree('tree1',"data: {0:d} channel, {1:d} samples".format(s1.nAdcCh, s1.nSamples))
    tree1.Branch('ch0',data1, "ch0[{0:d}][{1:d}]/F".format(s1.nAdcCh, s1.nSamples))

    for i in irange:
        s1.read_data([spattern.format(i),'xxx'], data1, data2)
        tree1.Fill()
    fout1.Write()

def compareWaveform(dsx):
    ### dsx = [(tag, file, chan, entry)]
    for a in dsx: print(a[1])
    fls = [TFile(a[1],'read') for a in dsx]
    trs = [f.Get('tree1') for f in fls]

    grs = [None]*len(dsx)
    opt = ''
    for i in range(len(dsx)):
#         print('adc[{0:d}]:Iteration$'.format(dsx[i][2]), "Entry$=={0:d}".format(dsx[i][3]), '')
#         print(trs[i].GetEntries())
#         print('adc[{0:d}]:Iteration$'.format(dsx[i][2]), "Entry$=={0:d}".format(dsx[i][3]), 'goff')
# #         trs[i].Draw('adc[{0:d}]:Iteration$'.format(dsx[i][2]), "Entry$=={0:d}".format(dsx[i][3]), 'goff')
# #         trs[i].Draw('adc[{0:d}]:Iteration$'.format(dsx[i][2]),"Entry$=={0:d}".format(dsx[i][3]), 'goff')
#         trs[i].Draw('adc[{0:d}]:Iteration$'.format(dsx[i][2]),"Entry$=={0:d}".format(dsx[i][3]), 'goff')
        print(i)
        trs[i].Draw('adc[{0:d}]:Iteration$'.format(dsx[i][2]),"Entry$=={0:d}".format(dsx[i][3]), opt)
        opt = 'same'
        print(i)
# #         grs[i] = gPad.GetPrimitive('Graph')
#         grs[i] = gDirectory.Get('Graph')
#         print(grs[i].GetN())
#         waitRootCmdX()
#         grs[i] = gDirectory.Get('htemp')
#         grs[i].SetName('gr'+str(i))
    
#     grs[0].Draw('AP')
#     for g in grs[1:]:
#         g.Draw('Psame PMC')
    waitRootCmdX()

def setPulse(v,f):
    cmd = 'ssh maydaq.dhcp.lbl.gov ./fungen_ctrl.py {0:.3f} {1:d}'.format(v,f)
    call(cmd, shell=True)

def take_calibration_samples(sTag, vs, n=5000):
    sc1 = SignalChecker()
    sc1.control_ip_port = "localhost:1024"
    dir1 = 'data/fpgaLin/'

    for v in vs:
#     for iv in range(22):
#         v = 0.025+iv*0.025
#         if v<0.175: continue
#         v = 0.05+iv*0.025
        setPulse(v,100)
        time.sleep(20)
        print "Taking sample with dU={0:.3f} mV".format(v*1000)
        sc1.take_samples2(n, dir1+sTag+"_{0:d}mV_f1000.root".format(int(v*1000)))

def take_data(sTag, n=5000, N=-1, dirx=None):
    sc1 = SignalChecker()
    sc1.control_ip_port = "localhost:1024"
    dir1 = 'data/fpgaLin/'

    ### put in a dedicated direcotry
    if dirx is not None:
        dir1 += dirx
        if not os.path.exists(dir1):
            os.makedirs(dir1)
    dir1 = dir1.rstrip()
    if dir1[-1] != '/': dir1 += '/'

    ### really start taking samples
    nSample = 0
    while nSample != N:
        print "Start sample {0:d}".format(nSample)
        status = sc1.take_samples2(n, dir1+sTag+"_data_{0:d}.root".format(nSample))
        if status: break
        nSample += 1

def test1():
    sc1 = SignalChecker()
    sc1.control_ip_port = "localhost:1024"
#     sc1.take_samples2(100, "data/test1.root")
#     sc1.take_samples2(5000, "data/sample1.root")
#     sc1.take_samples2(5000, "data/Jan18a_C2_50mV.root")
#     dir1 = 'data/fpgaLin/'
#     tag1 = dir1+'Jan22a_C2_100mV_'
#     for f in [100,200,500,1000]:
#         setPulse(0.1,f)
#         time.sleep(20)
#         sc1.take_samples2(5000, tag1+"f{0:d}.root".format(f))
#         sc1.check_enc2(tag1+"f{0:d}.root".format(f), tag1+"f{0:d}.dat".format(f))
    dir1 = 'data/fpgaLin/'
#     tag1 = dir1+'Jan22b_C2_'
#     for iv in range(16):
#         v = 0.025+iv*0.05
#         setPulse(v,1000)
#         time.sleep(20)
#         sc1.take_samples2(5000, tag1+"{0:d}mV_f1000.root".format(int(v*1000)))
#         sc1.check_enc2(tag1+"f{0:d}.root".format(f), tag1+"f{0:d}.dat".format(f))

#     tag1 = dir1+'Jan28b_C2_'
#     for iv in range(31):
#         v = 0.025+iv*0.025
# 
# #         if v<0.37: continue
#         setPulse(v,500)
#         time.sleep(50)
#         sc1.take_samples2(100, tag1+"{0:d}mV_f500.root".format(int(v*1000)))
#     sc1.take_samples2(200, dir1+"Jan31a_noise_dc.root")
#     sc1.take_samples2(5000, dir1+"Feb01a_noise_dc.root")
#     sc1.take_samples2(5000, dir1+"Feb05a_noise_dc_ch19.root")
    nSample = 0
    sTag = 'Feb25a'
#     while True:
#         sc1.take_samples2(1000, dir1+"Feb06c_data_{0:d}.root".format(nSample))
#         nSample += 1
    while True:
#         sc1.take_samples2(2000, dir1+"Feb08t1_data_{0:d}.root".format(nSample))
        print "Start sample {0:d}".format(nSample)
        sc1.take_samples2(1000, dir1+sTag+"_data_{0:d}.root".format(nSample))
        nSample += 1

#         if nSample == 1885: break
#     sc1.take_samples2(1000, dir1+"Feb06a_noise_dc_ch19.root")
#     sc1.take_samples2(5000, dir1+"Feb06b_noise_dc_ch19.root")

#
#     sc1.take_samples(10, name="Jan03a_{0:d}")
#     sc1.take_samples(5000, name="data/Jan04a/Jana04a_{0:d}")
#     sc1.take_samples2(5000, "data/Jan05a_150mV.root")
#     sc1.take_samples2(5000, "data/Jan05a_400mV.root")
#    sc1.take_samples2(5000, "data/Jan09a_300mV.root")
#     sc1.take_samples2(5000, "data/Jan05a_50mV.root")
#     sc1.take_samples2(5000, "data/Jan08a_100mV_R19p5312us.root")
#     sc1.take_samples2(5000, "data/Jan08a_100mV_R30p0us.root")
#     sc1.take_samples2(5000, "data/Jan08a_100mV_R40p0us.root")
#     sc1.take_samples2(5000, "data/Jan08a_100mV_r50p0us.root")
#     sc1.take_samples2(5000, "data/Jan08a_100mV_r30p0us.root")
#     sc1.take_samples2(5000, "data/Jan08a_100mV_r40p0us.root")
#     sc1.take_samples2(100, "data/test.root")
#     sc1.show_signal()
#     sc1.check_file('/data/Samples/TMSPlane/Dec26/sample_0.adc')
#     sc1.check_file('/data/Samples/TMSPlane/Dec27/Dec27a_1281.adc')
#     sc1.check_enc('/data/Samples/TMSPlane/Dec27/Dec27a_*.adc', ch=12)
#     sc1.check_enc2('data/root_files/Jan05a_50mV.root', 'Jan05a_50mV.dat')
#     sc1.check_enc2('data/root_files/Jan05a_100mV.root', 'Jan05a_100mV.dat')
#     sc1.check_enc2('data/root_files/Jan05a_150mV.root', 'Jan05a_150mV.dat')
#     sc1.check_enc2('data/root_files/Jan05a_400mV.root', 'Jan05a_400mV.dat')
#     sc1.check_enc2('data/root_files/Jan08a_100mV_r30p0us.root', 'Jan08a_100mV_r30p0us.dat')
#     sc1.check_enc2('data/root_files/Jan08a_100mV_r40p0us.root', 'Jan08a_100mV_r40p0us.dat')
#     sc1.check_enc2('data/root_files/Jan08a_100mV_r50p0us.root', 'Jan08a_100mV_r50p0us.dat')
#     sc1.check_enc2('/data/Samples/TMSPlane/root_files/Jan05a_50mV.root', 'Jan05a_50mV.dat')
#     sc1.check_enc2('/data/Samples/TMSPlane/root_files/Jan05a_150mV.root', 'Jan05a_150mV.dat')
#     sc1.check_enc2('data/sample1.root', 'Jan17_sample1.dat')
#     sc1.check_enc2('data/Jan18a_C2_50mV.root', 'data/Jan18a_C2_50mV.dat')
#     sc1.take_samples()
#     sc1.show_signal()
#     sc1.show_sample()
#     sc1.show_sample('/data/Samples/TMSPlane/Dec27/Dec27a_10{0:d}.adc',Ns=80,ich=12)
#     sc1.show_sample('/data/Samples/TMSPlane/Dec27/Dec27a_1000.adc',Ns=1,ich=12)

if __name__ == '__main__':
#     useLHCbStyle()
#     test1()
#     take_calibration_samples(sTag='Feb26a',n=5000)
#     take_calibration_samples(sTag='Mar07C1a',vs=[0.2+0.025*i for i in range(16)],n=3000)
#     take_calibration_samples(sTag='Feb26a',n=3000)
#     take_calibration_samples(sTag='Feb26b', vs=[0.025+0.05*i for i in range(10)],  n=3000)
#       take_data(sTag='Mar01t1b',n=200, N=1)
#       take_data(sTag='Mar04C1a',n=2000, N=1)
#       take_data(sTag='Mar05T1a',n=200, N=1)
#       take_data(sTag='Mar08D1a',n=1000, N=-1)
#       take_data(sTag='Apr04R1a',n=1000, N=-1)
#       take_data(sTag='Apr14T1a',n=1000, N=-1)
#       take_data(sTag='Apr22T1a',n=1000, N=-1)
      take_data(sTag='May10T1a',n=1000, N=-1, dirx='raw/May10T1a')
#     text2root(spattern='/data/Samples/TMSPlane/Dec27/Dec27a_{0:d}.adc',irange=range(10,20),outname='testxy.root')
#     text2root(spattern='data/Jan04a/Jana04a_{0:d}.adc',irange=range(5000),outname='ADC_Jan04a.root')
#     compareWaveform([('a1','/data/Samples/TMSPlane/root_files/Jan08a_100mV_R30p0us.root', 12, 20),('a2','/data/Samples/TMSPlane/root_files/Jan08a_100mV_R19p5312us.root', 12, 20)])
