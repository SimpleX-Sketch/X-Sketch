# X-Sketch

We define a new pattern in data streams: k-simplex items (k = 0, 1, 2), whose frequencies in several consecutive windows can be fitted by a polynomial of k-degree. We propose a novel sketch algorithm called X-Sketch for finding k-degree simplex items in real time, which is memory-efficient and accurate. Our experimental results show that the F1 Score of X-Sketch is on average 68.6%, 57.9%, and 42.2% higher than the baseline solution for k = 0, 1, 2, respectively. In addition, our case-specific experiments validate that it can also be applied to reduce the running time of the linear regression and ARIMA models mostly by more than 100 times.  
This paper has been submitted to IEEE ICDE 2023. We have uploaded the code of CPU and machine learning experiments (C++ version and Python version of X-Sketch), respectively, and put them in two folders.


### CPU

 You need to change the data format into the format of your data set in the code first

 cmake .

 make 

./cpu DATASET_PATH RESULT_PATH
