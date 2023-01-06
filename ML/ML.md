### ML

For running experiments on machine learning, you shall first enter directory cpu and modify Param.h. You shall modify line 5 to support predict mode. Then you can run cpu version of XSketch to get predict results of XSketch. It will also show items which can be predicted using machine learning. Finally, you need to enter directory ML and modify the variable `file` to input items to be predicted. Python will show the accuracy and predict time of two machine learning algorithms. 

### XSketch

For running Python experiments of XSketch, you need to compile dynamic libraries and change the data format in the code. The commands are as follows. Meanwhile, the dataset used by correct detector in cpp file shoulb be consistent with that used by python file.

```
g++ ../cpu/XSketch.cpp -o ./XSketch.so -g -Wall -O3 -std=c++17 -shared -fPIC

python3 XSketch.py
```

The accuracy and predict time of XSketch will be shown accordingly.
