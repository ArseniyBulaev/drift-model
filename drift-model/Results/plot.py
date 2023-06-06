import matplotlib.pyplot as plt
import numpy as np
import sys
import os

fname = sys.argv[1]
x, y = [], []
property_name = fname.split(".")[0]
print(os.getcwd())

with open(fname) as file:
    for line in file:
        x_value, y_value = line.split("\t")
        x.append(float(x_value))
        y.append(float(y_value))



fig, ax = plt.subplots()
ax.plot(x, y)

ax.set(xlabel='z', ylabel=property_name)
ax.grid()
plt.show()