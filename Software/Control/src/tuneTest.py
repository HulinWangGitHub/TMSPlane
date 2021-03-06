#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from PyDE import *
import threading
from Queue import Queue
import time
from TMS1mmX19Tuner import SensorConfig, CommonData
import socket
from command import *
import numpy as np
import matplotlib.pyplot as plt
import array
import json
from ROOT import *
gROOT.LoadMacro("sp.C+")
from ROOT import SignalProcessor
from reco_config import apply_config

class tuner(threading.Thread):
    def __init__(self, idx):
        threading.Thread.__init__(self)
        self.idx = idx
        self.rx_qs = None
        self.tx_qs = None
        self.atBounds = [(-10,10),(-10,10),(-10,10),(-10,10),(-10,10),(-10,10)] 
        self.atMaxIters = 1000

    def run(self):
        de = DE(self.auto_tune_fun, self.atBounds, maxiters=self.atMaxIters)
        ret = de.solve()
        self.tx_qs[self.idx] = None
        self.rx_qs[self.idx] = None
        print("Done {0:d}\n".format(self.idx))
 
        time.sleep(self.idx*1+1)
        print('-'*10, self.idx)
        print(ret)
       
    def auto_tune_fun1(self, x):
        return sum([(a-self.idx)*(a-self.idx) for a in x])

    def auto_tune_fun(self, x):
        ### put x to the queue
#         print(self.idx, 'sending', x)
        self.tx_qs[self.idx].put(x)
#         print('+++sent:['+','.join([str(a) for a in x]))

        ### wait for the returned value
        gg = self.rx_qs[self.idx].get()
#         print('+++get:{0:g}'.format(gg))
        return gg

class Train(threading.Thread):
    def __init__(self, cd):
        threading.Thread.__init__(self)
        self.cd = cd
        self.tx_qs = None
        self.rx_qs = None
        self.on = True
        self.mask = [0]*cd.nCh
        self.nSig = 3

        ### this copy of sensorVcodes is used to save the best value find so far
        self.sensorVcodes = [[v for v in cd.sensorVcodes[i]] for i in range(cd.nCh)]
#         self.bestConfigFile = 'current_best_config.json'
        self.bestConfigFile = 'C3_tt6_config.json'
        self.retBest = [0.]*cd.nAdcCh

        ### for data taking
        self.sc = None

        ### plotting
        self.axs = None
        self.plot_x = None
        self.pltCnt = 0
        self.pltN = 20

        self.sp = SignalProcessor()
        apply_config(self.sp, 'Helium')

        self.NVAL = 8
        s1 = self.cd.sigproc
        self.data1 = (s1.ANALYSIS_WAVEFORM_BASE_TYPE * (s1.nSamples * s1.nAdcCh))()
        self.ret1 = array.array('f',[0]*s1.nAdcCh)
        self.par1 = array.array('f',[0]*(self.cd.nCh*len(self.cd.inputVs)))
        self.t_values = array.array('f',[0]*(s1.nAdcCh*self.NVAL))
        self.keepAllData = False
        self.saveT0 = -1
        self.T = array.array('i',[0])
        self.Tag = array.array('i',[0])

    def setupOutput(self, outRootName='tt_test.root'):
        s1 = self.cd.sigproc
        self.fout1 = TFile(outRootName,'recreate')
        self.tree1 = TTree('tree1',"tune data for {0:d} channels, {1:d} samples".format(s1.nAdcCh, s1.nSamples))
        self.tree1.Branch('tag',self.Tag,'tag/I')
        self.tree1.Branch('T',self.T,'T/i')
        self.tree1.Branch('adc',self.data1, "adc[{0:d}][{1:d}]/F".format(s1.nAdcCh, s1.nSamples))
        self.tree1.Branch('ret',self.ret1, "ret[{0:d}]/F".format(s1.nAdcCh))
        self.tree1.Branch('par',self.par1, "par[{0:d}][{1:d}]/F".format(self.cd.nCh, len(self.cd.inputVs)))
        self.tree1.Branch('val',self.t_values, "val[{0:d}][{1:d}]/F".format(self.cd.nAdcCh, self.NVAL))

    def plot_data(self):
#         item = self.q.get()
        if self.axs is None:
            print("creating the plot......")
            import matplotlib
            matplotlib.rcParams['toolbar'] = 'None'
            plt.ion()
            fig, self.axs = plt.subplots(self.cd.nAdcCh, 1, sharex=True)
            # Remove horizontal space between axes
            plt.tight_layout(pad=0)
            fig.subplots_adjust(hspace=0)
            self.plot_x = [self.cd.adcDt * i for i in range(self.cd.nSamples)]

        for i in range(self.cd.nAdcCh):
            self.axs[i].cla()
#             self.axs[i].plot(self.cd.adcData[i])
#             self.axs[i].step(array.array('f', self.cd.adcData[i]), where='post')
#             self.axs[i].plot(self.plot_x, array.array('f', self.cd.adcData[i]))
            self.axs[i].plot(self.plot_x, self.data1[i*self.cd.nSamples : (i+1)*self.cd.nSamples])
        plt.draw()
#         self.q.task_done()


    def test_update_sensor(self, inputVs=None):
        print('/'*40)
        if inputVs is None:
            inputVs = [[1.379, 1.546, 1.626, 1.169, 1.357, 2.458],[1.379, 1.546, 1.626, 1.169, 1., 2.]]
        #### update sensor configurations
        ss = set()
        for i,p in enumerate(inputVs):
            if p is None: continue
            self.cd.set_sensor(i,p)
            ss.add(self.sc.tms1mmX19chainSensors[self.sc.tms1mmX19sensorInChain[i]][0])

        ### apply the configurations
        print(ss)
        for isr in ss:
            self.sc.update_sensor(isr)
        print('\\'*40)

    def take_data2(self, NEVT = 100):
        '''return an array of FOM'''
        s1 = self.cd.sigproc

        for ievt in range(NEVT):
            self.T[0] = int(time.time())
            self.cd.dataSocket.sendall(self.cd.cmd.send_pulse(1<<2));
            buf = self.cd.cmd.acquire_from_datafifo(self.cd.dataSocket, self.cd.nWords, self.cd.sampleBuf)
            s1.demux_fifodata(buf, self.data1, self.cd.sdmData)

            ### measure
            self.sp.measure_multipleX(self.data1, 2000, self.t_values)

            ### save
            self.tree1.Fill()

    def take_data(self):
        '''return an array of FOM'''
        s1 = self.cd.sigproc

        NV = 8
        NEVT = 100
        Values = [[0.]*(NV*NEVT) for i in range(self.cd.nAdcCh)]

        for ievt in range(NEVT):
            self.cd.dataSocket.sendall(self.cd.cmd.send_pulse(1<<2));
            buf = self.cd.cmd.acquire_from_datafifo(self.cd.dataSocket, self.cd.nWords, self.cd.sampleBuf)
            s1.demux_fifodata(buf, self.data1, self.cd.sdmData)

            self.sp.measure_multipleX(self.data1, 2000, self.t_values)
#             if ievt == 0:
#                 print("EVT",ievt) 
# #                 print([self.t_values[i] for i in range(NV*self.cd.nAdcCh)])
#                 print([self.t_values[6*NV+i] for i in range(NV)])

            for ich in range(self.cd.nAdcCh):
                for ipeak in range(NV):
                    Values[ich][ievt*NV+ipeak] = self.t_values[ich*NV+ipeak]

        for i in range(s1.nAdcCh):
#             self.ret1[i] = np.std(Values[i])/np.mean(Values[i])
            self.ret1[i] = -np.mean(Values[i])/np.std(Values[i])
            #             if i==6:
# #                 print(i,'->', Values[i][20*8:21*8])
#                 print(i,'->', Values[i][0*8:1*8])
#             print(i, np.std(Values[i]), np.mean(Values[i]), self.ret1[i])

        self.T[0] = int(time.time())
#         print(self.ret1[3],self.ret1[0])
#         for i in range(self.cd.nAdcCh): print(i,self.ret1[i])
        self.tree1.Fill()

        if self.T[0]-self.saveT0>200:
            self.tree1.AutoSave('SaveSelf')
            self.saveT0 = self.T[0]

    def run(self):
        nPar = len(self.cd.inputVs)

        ### save the initial values
        for i in range(self.cd.nCh):
            for kk in range(nPar): 
                self.par1[i*nPar+kk] = self.cd.tms1mmReg.dac_code2volt(self.sensorVcodes[i][kk])

        while self.on:
            cnt = 0
#             cnt1 = 0
            ss = set()
            for i,q in enumerate(self.rx_qs):
                if q is None: continue
                cnt += 1
#                 x = q.get()
#                 self.pars[i] = x

                if not q.empty():
#                     cnt1 += 1
                    x = q.get()
                    self.mask[i] = 1
#                     print('--- {0:d} get'.format(i)+' ['+','.join([str(a) for a in x])+']')
#                     self.pars[i] = x
#                     print(i,x)
                    self.cd.set_sensor(i,x)
                    for kk in range(nPar): self.par1[i*nPar+kk] = x[kk] 
                    ss.add(self.sc.tms1mmX19chainSensors[self.sc.tms1mmX19sensorInChain[i]][0])


            if cnt == 0: self.on = False
            if len(ss) == 0: continue
#             if cnt1 == 0: continue

            #### update sensor configurations
#             ss = set()
#             for i,p in enumerate(self.pars):
#                 if self.mask[i] == 1:
#                     self.cd.set_sensor(i,p)
#                     ss.append(self.sc.tms1mmX19chainSensors[self.sc.tms1mmX19sensorInChain[i]][0])

            ### apply the configurations
            for isr in ss:
                self.sc.update_sensor(isr,quiet=1)

            ### take data
            self.pltCnt += 1
            t = None
            if self.pltCnt%self.pltN == 0:
                self.plot_data()
#                 t = threading.Thread(target=self.plot_data)
#                 t.daemon = True
#                 t.start()
#                 self.q.put('run')

            time.sleep(30)
            if t is not None: t.join()

            self.take_data()
            ret = self.ret1

            ### return
            needUpdate = False
            for i,t in enumerate(self.tx_qs):
                if self.mask[i] == 1:
                    t.put(ret[i])
                    self.mask[i] = 0

                    ### save the values if it's the best so far
                    if ret[i]<self.retBest[i]:
                        print("find better parameters for channel", i)
                        print('old:',self.sensorVcodes[i], self.retBest[i])
                        self.retBest[i] = ret[i]
                        self.sensorVcodes[i] = [a for a in self.cd.sensorVcodes[i]]
                        print('new:',self.sensorVcodes[i], self.retBest[i])
                        needUpdate = True
#                     print("--- {0:d} {1:g}".format(i, self.meas[i]))
            if needUpdate:
                self.sc.write_config_fileX(self.bestConfigFile, self.sensorVcodes)
        print('Stopping the train.....')
        plt.close('all')
        self.tree1.Write()
        self.fout1.Close()

class TestClass:
    def __init__(self, nCh=19):
        self.x = None
        self.tx_qs = [None]*nCh
        self.rx_qs = [None]*nCh
        self.nCh = nCh
        self.muteList = []
        self.atBounds = None

    def prepare_train(self):
        host='192.168.2.3'
        if socket.gethostname() == 'FPGALin': host = 'localhost'

        dataIpPort = (host+':1024').split(':')
        sD = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        sD.connect((dataIpPort[0],int(dataIpPort[1])))

        ctrlIpPort = (host+':1025').split(':')
        sC = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        sC.connect((ctrlIpPort[0],int(ctrlIpPort[1])))

        cmd = Cmd()
        cd = CommonData(cmd, dataSocket=sD, ctrlSocket=sC)
        cd.aoutBuf = 1 # AOUT buffer select, 0:AOUT1, 1:AOUT2, >1:disable both
        cd.x2gain = 2 # BufferX2 gain 
        cd.sdmMode = 0 # SDM working mode, 0:disabled, 1:normal operation, 2:test with signal injection 
        cd.bufferTest = 0 #
#         cd.atTbounds = (2650,2750)
        cd.atTbounds = (4050,4150)

        sc1 = SensorConfig(cd, configFName='C3_config.json')

        tr1 = Train(cd)
        tr1.tx_qs = self.rx_qs
        tr1.rx_qs = self.tx_qs
        tr1.sc = sc1

        return tr1

    def save_config(self, cList, oName, fcName = 'Mar05Te_tt_test.root'):
        ### base on the default configuration, overwite it with new ones
        inputVs = [None]*self.nCh

        ### get files
        fin = TFile(fcName,'read')
        ch = fin.Get('tree1')

        ## get the parameters
        cd = CommonData(cmd=None)
        for ich in range(self.nCh):
            if cList[ich] is None: continue

            ievt = cList[ich]
            n1 = ch.Draw('par[{0:d}]:ret[{0:d}]'.format(ich),'Entry$=={0:d}'.format(ievt),'goff')
            v1 = ch.GetV1()
            cd.set_sensor(ich,[v1[j] for j in range(n1)])

        ### save
        config = {}
        for i in range(cd.nCh):
            config[i] = dict(zip(cd.voltsNames, cd.sensorVcodes[i]))
        with open(oName, 'w') as fp:
            fp.write(json.dumps(config, sort_keys=True, indent=4))

    def save_config1(self, cList0, oName, fcName = 'Mar05Te_tt_test.root'):
        ### base on the default configuration, overwite it with new ones
        inputVs = [None]*self.nCh

        ### get files
        fin = TFile(fcName,'read')
        ch = fin.Get('tree1')

        cut = ''
        ### get the entry numbers
        cList = [None]*self.nCh
        for ich in range(self.nCh):
            if cList0[ich] is None: continue

            n = ch.Draw("ret[{0:d}]:Entry$".format(ich),cut,"goff")
            v1 = ch.GetV1()
            v2 = ch.GetV2()
            vx = sorted([(v1[i],int(v2[i])) for i in range(n)], key=lambda x:x[0])
            cList[ich] = vx[cList0[ich]][1]

        print(cList)
        self.save_config(cList, oName, fcName)


    def validate_tune(self, fcName='C3_tt3.root', oName='C3_tt3_valid0.root'):
        '''Here will will use low frequency pulse'''
#         fcName = 'Mar05Te_tt_test.root'
#         oName = 'Mar05Te_tt_valid4.root'
#         fcName = 'C3_tt3.root'
#         oName = 'C3_tt3_valid0.root'
        dT_wait = 50
        N_data = 1000
        topN = 10

        cList = [None]*self.nCh
        ### get the best config
        fin = TFile(fcName,'read')
        ch = fin.Get('tree1')

        cut = ''
        for ich in range(self.nCh):
            n = ch.Draw("ret[{0:d}]:Entry$".format(ich),cut,"goff")
            v1 = ch.GetV1()
            v2 = ch.GetV2()
            vx = sorted([(v1[i],int(v2[i])) for i in range(n)], key=lambda x:x[0])
            cList[ich] = vx[:topN]

        tr1 = self.prepare_train()
        tr1.setupOutput(oName)
        nPar = len(tr1.cd.inputVs)

        ### take the default one first
        tr1.Tag[0] = -1
        time.sleep(dT_wait)
        for ich in range(self.nCh):
            for j in range(nPar):
                tr1.par1[ich*nPar+j] = tr1.cd.tms1mmReg.dac_code2volt(tr1.cd.sensorVcodes[ich][j])
        tr1.take_data2(N_data)

        ### check those from tune
        for topi in range(topN):
            inputVs = [None]*self.nCh
            ## get the parameters
            for ich in range(self.nCh):
                ievt = cList[ich][topi][1]
                n1 = ch.Draw('par[{0:d}]:ret[{0:d}]'.format(ich),'Entry$=={0:d}'.format(ievt),'goff')
                v1 = ch.GetV1()
                inputVs[ich] = [v1[j] for j in range(n1)]
                for j in range(nPar): tr1.par1[ich*nPar+j] = v1[j]
                tr1.ret1[ich] = ch.GetV2()[0]

            ## take data
            tr1.test_update_sensor(inputVs)
            tr1.Tag[0] = topi
            time.sleep(dT_wait)
            tr1.take_data2(N_data)

        ## save the tree and close
        tr1.tree1.Write()
        tr1.fout1.Close()

    def test_tune(self, oName='tt_out1.root'):
        tr1 = self.prepare_train()
        tr1.on = True
        tr1.setupOutput(oName)

#         tr1.test_update_sensor()
        tr1.take_data()
#         return
        ### for the chip #3? the default one in LBL
#         better_bounds = [None]*self.nCh
#         better_bounds[0] = [(0.9, 1.5), (0.5, 1.2), (1.1, 1.8), (0.5, 1.4), (1.4, 1.8), (2.5, 2.8)]
#         better_bounds[1] = [(0.8, 1.5), (0.5, 1.4), (1.0, 1.8), (0.5, 1.6), (1.4, 2.0), (2.5, 2.8)]
#         better_bounds[2] = [(0.5, 1.5), (0.5, 1.8), (0.5, 1.5), (0.5, 2.0), (0.8, 2.4), (1.8, 2.8)]
#         better_bounds[3] = [(0.8, 1.5), (0.5, 1.3), (1.0, 1.8), (0.5, 2.0), (0.8, 2.0), (1.8, 2.8)]
#         better_bounds[4] = [(0.9, 1.2), (0.5, 1.3), (0.5, 1.8), (0.5, 1.4), (0.8, 1.8), (2.5, 2.8)]
#         better_bounds[5] = [(0.8, 1.5), (0.5, 1.5), (0.9, 1.8), (0.5, 2.0), (0.8, 2.0), (1.8, 2.8)]
#         better_bounds[6] = [(0.8, 1.2), (0.5, 1.5), (1.0, 1.5), (0.5, 2.0), (0.8, 1.8), (1.8, 2.8)]
#         better_bounds[7] = [(0.5, 1.5), (0.5, 1.1), (1.4, 1.8), (0.5, 1.6), (0.8, 1.5), (2.5, 2.8)]
#         better_bounds[8] = [(0.5, 1.5), (0.5, 1.8), (0.5, 1.8), (0.5, 2.0), (0.8, 1.6), (2.4, 2.8)]
#         better_bounds[9] = [(0.7, 1.5), (0.5, 1.4), (1.0, 1.8), (0.5, 1.8), (0.8, 2.0), (2.5, 2.8)]
#         better_bounds[10] = [(0.8, 1.3), (0.5, 1.4), (0.9, 1.6), (0.5, 2.0), (0.8, 2.0), (1.9, 2.8)]
#         better_bounds[11] = [(0.5, 1.5), (0.5, 1.0), (1.3, 1.8), (0.5, 1.4), (0.8, 1.2), (2.5, 2.8)]
#         better_bounds[12] = [(0.7, 1.3), (0.5, 1.3), (1.1, 1.8), (0.5, 1.8), (0.8, 2.0), (2.5, 2.8)]
#         better_bounds[13] = [(0.8, 1.3), (0.5, 1.4), (1.1, 1.8), (0.5, 1.7), (0.8, 2.2), (2.4, 2.8)]
#         better_bounds[14] = [(0.7, 1.5), (0.5, 1.2), (1.3, 1.8), (0.5, 1.5), (0.8, 1.8), (2.5, 2.8)]
#         better_bounds[15] = [(0.9, 1.3), (0.5, 1.2), (1.1, 1.8), (0.5, 1.6), (0.8, 1.8), (2.5, 2.8)]
#         better_bounds[16] = [(0.9, 1.5), (0.5, 1.2), (1.2, 1.8), (0.5, 1.5), (0.8, 1.8), (2.5, 2.8)]
#         better_bounds[17] = [(0.5, 1.5), (0.5, 1.0), (1.0, 1.5), (0.5, 1.1), (0.8, 1.2), (2.5, 2.8)]
#         better_bounds[18] = [(0.9, 1.4), (0.5, 1.3), (1.0, 1.8), (0.5, 1.6), (0.8, 1.8), (2.5, 2.8)]

        ### for the chip #3? the default one in LBL
        better_bounds = [None]*self.nCh
        better_bounds[0]  = [(1.3, 1.6), (0.7, 1.4), (1.5, 2.0), (0.6, 1.7), (1.2, 1.5), (2.2, 2.5)]  # 0 
        better_bounds[1]  = [(1.0, 1.7), (0.6, 1.3), (0.8, 1.5), (0.6, 1.7), (0.9, 1.7), (2.1, 2.5)]  # 1
        better_bounds[2]  = [(1.0, 1.7), (0.6, 1.5), (0.9, 1.7), (0.6, 1.7), (1.1, 1.8), (2.3, 2.6)]  # 2
        better_bounds[3]  = [(1.3, 1.7), (0.6, 1.6), (1.5, 2.0), (0.6, 1.7), (1.0, 1.5), (2.1, 2.5)]  # 3
#         self.atBounds = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (0.9, 1.7), (2.3, 2.6)]
        better_bounds[4]  = [(1.0, 1.5), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (0.9, 1.7), (2.2, 2.6)]  # 4
        better_bounds[5]  = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.5), (0.6, 1.7), (0.8, 1.6), (2.1, 2.5)]  # 5
        better_bounds[6]  = [(1.2, 1.7), (0.6, 1.6), (1.5, 2.0), (0.6, 1.7), (0.9, 1.6), (2.2, 2.5)]  # 6
        better_bounds[7]  = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (0.9, 1.7), (2.3, 2.5)]  # 7
        better_bounds[8]  = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.5), (0.6, 1.7), (0.9, 1.7), (2.1, 2.5)]  # 8
        better_bounds[9]  = [(1.0, 1.6), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (1.1, 1.9), (2.3, 2.8)]  # 9
        better_bounds[10] = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (0.9, 1.7), (2.1, 2.8)]  #10
        better_bounds[11] = [(1.1, 1.6), (0.6, 1.6), (1.5, 2.0), (0.6, 1.7), (1.0, 1.5), (2.1, 2.5)]  #11 
        better_bounds[12] = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.5), (0.6, 1.7), (0.9, 1.6), (2.3, 2.6)]  #12 
        better_bounds[13] = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (0.9, 1.7), (2.1, 2.8)]  #13 
        better_bounds[14] = [(1.0, 1.7), (0.6, 1.6), (0.9, 1.7), (0.6, 1.7), (0.5, 1.5), (2.4, 2.8)]  #14 
        better_bounds[15] = [(1.1, 1.7), (0.6, 1.6), (0.9, 1.8), (0.6, 1.7), (0.9, 1.7), (2.3, 2.6)]  #15 
        better_bounds[16] = [(1.0, 1.7), (0.6, 1.6), (0.5, 1.7), (0.6, 1.7), (0.9, 1.7), (2.1, 2.8)]  #16 
        better_bounds[17] = [(1.3, 1.6), (0.6, 1.6), (0.5, 2.0), (0.6, 1.7), (1.2, 1.6), (2.1, 2.5)]  #17 
        better_bounds[18] = [(1.2, 1.6), (0.6, 1.6), (1.5, 2.0), (0.6, 1.7), (1.0, 1.6), (2.1, 2.5)]  #18 

        for i in range(self.nCh):
            if i in self.muteList: continue
            self.tx_qs[i] = Queue()
            self.rx_qs[i] = Queue()
            th1 = tuner(i)
            th1.tx_qs = self.tx_qs
            th1.rx_qs = self.rx_qs            
#             th1.atBounds = tr1.cd.atBounds
            th1.atBounds = better_bounds[i]
#             th1.atBounds = better_bounds[i] if better_bounds[i] is not None else cd.atBounds
            th1.atMaxIters = tr1.cd.atMaxIters
            th1.start()
        tr1.start()

def test1():
    tc1 = TestClass()
#     tc1.muteList = [3,5,6,8,12,18]
#     tc1.muteList = []
#     tc1.test_tune('C3_tt7.root')
#     tc1.validate_tune(fcName='C3_tt6.root', oName='C3_tt6_valid1.root')
#     cList0 = [   8,   2,    2,    7,    3,   0,    0,    0,   2,    1,    2,    0,    1,    1,    0,   2,    4,   2,     2]
#     cList = [2148, 650, 1288, 2294, 1580, 420, 1521, 1297, 509, 1100, 1258, 1186, 1762, 1703, 1747, 513, 1750, 790, 2192 ]
#     cList0 = [   8,   4,    2,    7,    5,   0,    0,    0,   2,    1,    2,    0,    1,    3,    0,   2,    4,    1,    2]
#     cList = [2148, 670, 1288, 2294, 2299, 420, 1521, 1297, 509, 1100, 1258, 1186, 1762,  901, 1747, 513, 1750, 1842, 2192 ]
#     cList0= [   8,   4,    2,    7,    5,   0,    0,    0,   2,    3,    2,    0,    1,    3,    0,   2,    4,    x,    2]
#     cList = [2148, 670, 1288, 2294, 2299, 420, 1521, 1297, 509, 1870, 1258, 1186, 1762,  901, 1747, 513, 1750, 1343, 2192 ]
#     cList = [2148, 670, 804, 2294, 2299, 420, 1521, 1297, 509, 1870, 1258, 1186, 1762,  901, 1747, 513, 1750, 1343, 2192 ]
#     tc1.save_config(cList,'new_config.json',fcName='Mar05Te_tt_test.root')
#     tc1.save_config1([0]*tc1.nCh,'new_C3_config.json',fcName='C3_tt3.root')
    elist = [0]*tc1.nCh
#     elist[4] = 
    elist[6] = 6

#     elist[7] = 
#     elist[10] = 
    elist[11] = 2

#     elist[13] = 
#     elist[14] = 
    elist[15] = 2
#     elist[16] = 
    tc1.save_config1(elist,'new_C3_config6.json',fcName='C3_tt6.root')

if __name__ == '__main__':
    test1()
