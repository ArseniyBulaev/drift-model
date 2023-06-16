import matplotlib.pyplot as plt
import numpy as np
import sys
import os

fnames = sys.argv[1:]
X, Y = [], []
print(fnames)
print(X)
property_name = fnames[0].split(".")[0]

for i, fname in enumerate(fnames):
    X.append([])
    Y.append([])
    with open(fname) as file:
        for line in file:
            x_value, y_value = line.split("\t")
            X[i].append(float(x_value))
            Y[i].append(float(y_value))



fig, ax = plt.subplots()
for x, y in zip(X, Y):
    ax.plot(x, y)
    
ax.set(xlabel='z', ylabel=property_name)
ax.grid()
plt.show()