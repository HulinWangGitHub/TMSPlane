#!/usr/bin/env python3
from multiprocessing import Pool
import matplotlib.pyplot as plt
# from scipy import signal
from scipy.signal import find_peaks, peak_prominences
import numpy as np
from sys import argv
from glob import glob
import os, time, re
import ROOT
from root_numpy import root2array, array2root

#signIt = lambda t: t if t<127 else t-256
signIt = lambda t: t 

myOutPutDir = './'

class findrate(object):
    def __init__(self, filename="tek0001CH2.isf"):
        self.filename = filename
        self.timelap = 10. #second
        self.inputarray = np.array([])
        self.npeaks = 0
        self.date = ""
        self.time = ""
        self.diffpeakmean = 0
        self.heightmean = 0
        self.prominencemean = 0

    def getinput(self):
        with open(self.filename,'rb') as f1:
            b0 = f1.read(1024)
            b1i = b0.find(b':CURVE')
            b = b0[:b1i]
            c = b.split(b';')

            b2 = b0[b1i:]
            c += [b2]

            pars = {}
            for v in c:
                a = v.find(b" ")
                print(v,a)
                pars[v[:a]] = v[a+1:]
    #         print(pars)

#             Nbyte = pars[b'BIT_NR']
#             print(type(Nbyte),Nbyte)
#             return

            self.timelap = float(pars[b'WFID'].split(b',')[3].split(b'/')[0][1:6])
            self.timelap *=10
            self.date = pars[b':DATE'].decode("utf-8")
            self.time = pars[b':TIME'].decode("utf-8")
    
            wav = pars[b':CURVE'][10:-1]
            wav += f1.read()
            wav1 = [signIt(int(x)) for x in wav[::2]] if pars[b'BIT_NR'] == b'16' else [signIt(int(x)) for x in wav]
    #         wav1 = [from_bytes(x,'big',signed=True) for x in wav[::2]]
    #        wav1 = [ax(j) for j in wav[::2]]
            print(len(wav1), wav1[:10])

            self.inputarray = np.asarray(wav1)

    def processinput(self, N):
        basename = os.path.basename(self.filename)
        arwav2 = self.inputarray[:N] if N>0 else self.inputarray[:]
        self.inputarray = None
#         peaks = signal.find_peaks(wav2) 
        #peaks, _ = find_peaks(arwav2, height=8, width=150, distance=850, prominence=None) 
        peaks, properties = find_peaks(arwav2, height=None, width=100, prominence=10, distance=500) 

        print(peaks)
        #wav3 = [wav2[ix] for ix in peaks]

        promsall = peak_prominences(arwav2, peaks, wlen=None)
        proms = promsall[0]
        promxmins = promsall[1]
        promxmaxs = promsall[2]
        lpromscallow = []
        for i in range(len(proms)-1):
            if arwav2[promxmins[i]] < arwav2[promxmaxs[i]]:
                lpromscallow.append(arwav2[peaks[i]] - arwav2[promxmins[i]])
            else:
                lpromscallow.append(arwav2[peaks[i]] - arwav2[promxmaxs[i]])
        promscallow = np.asarray(lpromscallow)
        print(proms[:10])
        print(promscallow[:10])



        contour_heights = arwav2[peaks] - proms
        plt.plot(arwav2)
        plt.plot(peaks, arwav2[peaks], "x")
        plt.vlines(x=peaks, ymin=contour_heights, ymax=arwav2[peaks], color = "C1")
        plt.vlines(promxmins, ymin=contour_heights, ymax=arwav2[peaks], color = "C2", linestyles = "dashed")
        plt.vlines(promxmaxs, ymin=contour_heights, ymax=arwav2[peaks], color = "C2", linestyles = "dotted")
        plt.hlines(y=arwav2[peaks], xmin=promxmins, xmax=promxmaxs, color = "C2", linestyles = "dashdot")
        plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"], xmax=properties["right_ips"], color = "C1")

        #x2 = [peaks[i+1]-peaks[i] for i in range(len(peaks)-1)]
        fig0, ax0 = plt.subplots(num="diffpeak_"+basename)

        num_bins = 50
        # the histogram of the data
        n, bins, patches = ax0.hist(np.diff(peaks), num_bins)
        

        fig1, ax1 = plt.subplots(num="peakheight_"+basename)
        #n, bins, patches = ax.hist([arwav1[ix] for ix in peaks], num_bins)
        ax1.hist(arwav2[peaks], num_bins)

        fig2, ax2 = plt.subplots(num="peakprom_"+basename)
        #plt.plot(proms)
        ax2.hist(proms, num_bins)
        #plt.show()

        self.npeaks = len(peaks)
        self.diffpeakmean = np.mean(np.diff(peaks))
        self.heightmean = np.mean(arwav2[peaks])
        self.prominencemean = np.mean(proms)
#         artmp = np.column_stack((proms, arwav2[peaks]))
#         print(dir(artmp))
#         artmp.dtype = [['a', np.float32], ['b', np.float32]]
        #rootfile = array2root(np.array(zip(proms, arwav2[peaks]), dtype=[('a', np.float32), ('b', np.float32)])  , basename + '.root',  mode="recreate")
        #rootfile = array2root( artmp, basename + '.root',  mode="recreate")
#         rootfile = array2root( np.array(np.column_stack((proms, arwav2[peaks])), dtype=[('a', np.float32), ('b', np.float32)])  , basename + '.root',  mode="recreate")
        df = np.append(np.array([-1]),np.diff(peaks))
        pks = arwav2[peaks]
        bigx = np.array([(df[i],pks[i],proms[i]) for i in range(len(df))],dtype=[('dt',np.float32),('pk',np.float32),('proms',np.float32)])
        rootfile = array2root( bigx, myOutPutDir+"/"+basename[:-4] + '.root',  mode="recreate")

def processor(inputName):
    myrate = findrate(inputName)
    myrate.getinput()
    myrate.processinput(-1)
    return myrate

def multi_run(mydir='./test1'):
    script, pattern, Nfiles = argv

    if not os.path.exists(mydir): os.makedirs(mydir)
    global myOutPutDir
    myOutPutDir = mydir

    filesall = sorted([f for f in glob(pattern) if not os.path.exists(mydir+'/'+os.path.basename(f)[:-4]+'.root')], key=lambda f:os.path.getmtime(f))
    files = filesall[:int(Nfiles)] if int(Nfiles)>0 else filesall[:]
    print(files)
    myrates = []
    print(len(files))

    p = Pool(6)
    myrates = p.map(processor, files)
   
    newfile = open(mydir+"/summary.txt","a")
    print("filename", "date", "time", "npeaks", "duration", "rates", "heightmean", "prominencemean")
    #newfile.write(" ".join(["filename", "date", "time", "npeaks", "duration", "rates", "heightmean", "prominencemean"]))
    for i in range(len(myrates)):
        print(os.path.basename(myrates[i].filename), myrates[i].date, myrates[i].time, myrates[i].npeaks, myrates[i].timelap, myrates[i].npeaks/myrates[i].timelap, myrates[i].heightmean, myrates[i].prominencemean)
        newfile.write(f"{os.path.basename(myrates[i].filename)} {myrates[i].date} {myrates[i].time} {myrates[i].npeaks} {myrates[i].timelap} {myrates[i].npeaks/myrates[i].timelap} {myrates[i].heightmean} {myrates[i].prominencemean}\n")   


def test():
    script, pattern, Nfiles = argv
    filesall = sorted([f for f in glob(pattern)], key=lambda f:os.path.getmtime(f))
    files = filesall[:int(Nfiles)] if int(Nfiles)>0 else filesall[:]
    print(files)
    myrates = []
    print(len(files))
    for i in range(len(files)):
        print(f"processing the {i} of {len(files)} file: {files[i]}")
        myrate = findrate(files[i])
        myrate.getinput()
        myrate.processinput(-1)
        myrates.append(myrate)
    print("filename", "date", "time", "npeaks", "duration", "rates", "heightmean", "prominencemean")
    for i in range(len(myrates)):
        print(os.path.basename(myrates[i].filename), myrates[i].date, myrates[i].time, myrates[i].npeaks, myrates[i].timelap, myrates[i].npeaks/myrates[i].timelap, myrates[i].heightmean, myrates[i].prominencemean)
    plt.show()

    #summary plots
    #ROOT.TH1F 

if __name__ == '__main__':
    test()
    #multi_run()
