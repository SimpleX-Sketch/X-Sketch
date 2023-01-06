import numpy as np
import pandas as pd
import math
import os
import time
from tqdm import tqdm
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import adfuller as ADF  
from statsmodels.graphics.tsaplots import plot_pacf 
from statsmodels.stats.diagnostic import acorr_ljungbox 
from statsmodels.tsa.arima.model import ARIMA
import warnings
warnings.filterwarnings("ignore")

# parameters for linear regressing
window_size = 10000
run_length = 300000
last_timestamp = 0
win_cnt = 0
record = {}
dic = {}
predict = {}

# test
correct = 0
total = 0

# file input and get simplex items
path = "/share/datasets/CAIDA2016/1.dat"
file = "result.txt"
f = open(path, 'rb')
g = open(file)
t = int(g.readline())
for _ in range(t):
	i, j = map(int, g.readline().split())
	# print(i, j)
	if i in dic:
		dic[i].append(j)
	else:
		dic[i] = [j] 

# read the data
ids = []
for _ in range(run_length):
	a = int.from_bytes(f.read(8), byteorder='little', signed=False)
	b = int.from_bytes(f.read(8), byteorder='little', signed=False)
	ids.append(b)

# while len(ids) < run_length:
# 	l = f.readline().split()
# 	for i in l:
# 		ids.append(int(i))
f.close()
g.close()

start = time.time()
for index in tqdm(range(run_length)):
	if index >= last_timestamp + window_size:
		if (win_cnt + 1) in dic and win_cnt > 6:
			predict_result = {}
			for i in dic[win_cnt + 1]:
				# using time-series to predict
				bic_matrix = []
				arr = DataFrame(record[i])
				for p in range(3):
					tmp = []
					for q in range(3):
						try:
							tmp.append(ARIMA(arr, order = (p,1,q)).fit().bic)
						except:
							tmp.append(None)
					bic_matrix.append(tmp)
				bic_matrix = pd.DataFrame(bic_matrix)
				# print(bic_matrix)
				try:
					p, q = bic_matrix.stack().idxmin()
				except:
					continue
				model = ARIMA(arr, order = (p,1,q)).fit()
				predict_result[i] = model.forecast(steps=1)[win_cnt + 1]
			predict[win_cnt + 1] = predict_result
		if win_cnt in dic and win_cnt > 7:
			for i in dic[win_cnt]:
				if len(record[i]) <= win_cnt:
					truth = 0
				else:
					truth = record[i][win_cnt]
				if i not in predict[win_cnt]:
					pr = 0
				else:
					pr = predict[win_cnt][i]
				if abs(truth - pr) <= 5 or (truth * 2 >= pr and truth <= 2 * pr):
					correct += 1
				total += 1
		win_cnt += 1
		last_timestamp += window_size
	if ids[index] not in record:
		record[ids[index]] = [0] * (win_cnt + 1)
	while len(record[ids[index]]) <= win_cnt:
		record[ids[index]].append(0)
	record[ids[index]][win_cnt] += 1
end = time.time()

print("correct: %d total: %d\naccuracy: %.4f" % (correct, total, 100.0 * correct / total))
print("time: %.2f" % (end - start))
print("throughput: %.2f" % (run_length / (end - start) / 1000000))