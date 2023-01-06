import os
from ctypes import *
import math
import time
import numpy as np
from time import perf_counter

lib = CDLL('./XSketch.so')
lib.XSketch_new.restype = c_void_p
lib.XSketch_new.argtypes = [c_int32, c_double, c_double, c_double, c_int32, c_int32, c_double]
lib.insert.restype = None
lib.insert.argtypes = [c_void_p, c_ulonglong, c_int32]
lib.query.restype = py_object
lib.query.argtypes = [c_void_p]
lib.report.restype = py_object
lib.report.argtypes = [c_void_p]
lib.eval.restype = None
lib.eval.argtypes = [c_void_p, c_double, c_double, c_int32]

# sm.query()
# [id1, f1, id2, f2, id3, f3, ...]
# sm.report()
# [id1, start1, end1, id2, start2, end2, ...]

class XSketch:
    def __init__(self, lib, memory, var_input, error_input, ratio_input, CellNum, S1Wins, potential):
        self.lib = lib
        self.memory = memory
        self.var_input = var_input
        self.error_input = error_input
        self.ratio_input = ratio_input
        self.CellNum = CellNum
        self.S1Wins = S1Wins
        self.potential = potential
        self.obj = self.lib.XSketch_new(memory, var_input, error_input, ratio_input, CellNum, S1Wins, potential)
    
    def insert(self, id, timestamp):
        self.lib.insert(self.obj, id, timestamp)

    def query(self):
        return self.lib.query(self.obj)
    
    def report(self):
        return self.lib.report(self.obj)
    
    def eval(self, run_length):
        self.lib.eval(self.obj, self.var_input, self.error_input, run_length)



if __name__ == '__main__':

    run_length = 300000
    ids = []
    
    # CAIDA16
    f = open('/share/datasets/CAIDA2016/1.dat', 'rb')
    for _ in range(run_length):
        a = int.from_bytes(f.read(8), byteorder='little', signed=False)
        b = int.from_bytes(f.read(8), byteorder='little', signed=False)
        ids.append(b)
    
    # Transactional
    # with open('/share/gjr/T10I4D100K.dat','r') as f:
    #     lines = f.readlines()
    #     i = 0
    #     while i <= run_length:
    #         for line in lines:
    #             line = line.split(' ')[:-1]
    #             for a in line:
    #                 ids.append(int(a))
    #                 i += 1
    #                 if i == run_length:
    #                     break
    #         if i == run_length:
    #             break
    
    # MAWI
    # f = open('/share/lxd/time07.dat', 'rb')
    # for _ in range(run_length):
    #     a = int.from_bytes(f.read(8), byteorder='little', signed=False)
    #     b = int.from_bytes(f.read(16), byteorder='little', signed=False)
    #     ids.append(a)
    
    # Data Center/ Synthetic
    # f = open('/share/lxd/dc.dat', 'rb')
    # f = open('/share/zipf_2022/zipf_1.5.dat', 'rb')
    # for _ in range(run_length):
    #     a = int.from_bytes(f.read(4), byteorder='little', signed=False)
    #     ids.append(a)
    
    sm = XSketch(lib, 300, 0.1, 0.5, 0.8, 4, 3, 0.25)
    
    start = perf_counter()
    for i in range(run_length):
        sm.insert(ids[i], i)
    finish = perf_counter()
    time.sleep(1)

    sm.eval(run_length)
    print('Run time: {}s'.format(finish - start))
    # print('Throughput: ', run_length / (1.0 * (finish - start) * 1000))

