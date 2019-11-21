#!/usr/bin/env python
from filterChecker import filterChecker

fc1 = filterChecker()
fc1.filterPars = [100,200,1000,-1]
fc1.recoCfg = 'Lithium'
fc1.chan = 8
# fc1.offline_check('/data/Samples/TMSPlane/fpgaLin/raw/Nov04c/Nov04c_100mV_data_0.root')
fc1.offline_check('/data/Samples/TMSPlane/fpgaLin/raw/Nov04c/Nov04c_40mV_data_0.root')
