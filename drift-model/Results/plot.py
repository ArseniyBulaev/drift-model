import matplotlib.pyplot as plt
import numpy as np
import sys
import os

fnames = [line.rstrip() for line in sys.stdin]
print(fnames)
X, Y = [], []
property_name = fnames[0].split("__")[0]
times = [fname.split("__")[1].split(".")[0] for fname in fnames]


for i, fname in enumerate(fnames):
    X.append([])
    Y.append([])
    with open(fname) as file:
        for line in file:
            x_value, y_value = line.split("\t")
            X[i].append(float(x_value))
            Y[i].append(float(y_value))



fig, ax = plt.subplots()
for x, y, time in zip(X, Y, times):
    ax.plot(x, y, label = time + " sec.")
 
ax.legend() 
ax.set(xlabel='z', ylabel=property_name)
ax.grid()
plt.show()