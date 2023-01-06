import numpy as np
import pandas as pd
import math
import os
import time
from tqdm import tqdm

# parameters for linear regressing
P = 5
K = 2
window_size = 10000
run_length = 300000
last_timestamp = 0
win_cnt = 0
record = {}
dic = {}
predict = {}
x = np.array(range(P))

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
f.close()
g.close()

start = time.time()
for index in tqdm(range(run_length)):
	if index >= last_timestamp + window_size:
		# window ends, start linear regressing
		if win_cnt >= P - 1:
			predict_result = {}
			for i in record:
				y = []
				for j in range(P):
					y.append(record[i][(win_cnt + 1 + j) % P])
				y = np.array(y)
				beta = np.polyfit(x, y, K)
				poly = np.poly1d(beta)
				p = poly(P)
				predict_result[i] = max(0, round(p))
			predict[win_cnt + 1] = predict_result
			# start possible testing
			if win_cnt in dic:
				for i in dic[win_cnt]:
					machine_learning = predict[win_cnt][i]
					ground_truth = record[i][win_cnt % P]
					total += 1
					if abs(machine_learning - ground_truth) <= 5 or (machine_learning * 2 >= ground_truth and machine_learning < 2 * ground_truth):
						correct += 1
			# clear the window
			for i in list(record.keys()):
				record[i][(win_cnt + 1) % P] = 0
				if record[i] == [0] * P:
					record.pop(i)
		win_cnt += 1
		last_timestamp += window_size
	if ids[index] not in record:
		record[ids[index]] = [0] * P
	record[ids[index]][win_cnt % P] += 1
end = time.time()

print("correct: %d total: %d\naccuracy: %.4f" % (correct, total, 100.0 * correct / total))
print("time: %.2f" % (end - start))
print("throughput: %.2f" % (run_length / (end - start) / 1000000))