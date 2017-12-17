import numpy as np
import os
import matplotlib.pyplot as plt


fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Energy graph")
ax1.set_xlabel('Iterations')
ax1.set_ylabel('Energy')
x = []
lre = []
diff = []
abs_diff = []
first = True
with open("Energy_graph.txt") as f:
    lines = f.readlines()
    for line in lines:
        if line==lines[0]:
            if len(x)>40 and first:
                first = False
                ax1.plot(x, lre, c='r', label='energy curve')
                ax1.plot(x, diff, c='b', label='difference energy curve')
                ax1.plot(x, abs_diff, c='g', label='absolute difference')
            if len(x)>40 and not first:
                ax1.plot(x, lre, c='r')
                ax1.plot(x, diff, c='b')
                ax1.plot(x, abs_diff, c='g')
            x = []
            lre = []
            diff, abs_diff = [], []
        else : 
            try:
                l = line.split(";")
       
                x.append(float(l[0]))
                lre.append(float(l[1]))
                diff.append(float(l[2]))
                abs_diff.append(abs(float(l[2])))
            except:
               print(line)

if len(x)>20 and first:
    first = False
    ax1.plot(x, lre, c='r', label='energy curve')
    ax1.plot(x, diff, c='b', label='difference energy curve')
    ax1.plot(x, abs_diff, c='g', label='absolute difference')
if len(x)>40 and not first:
    ax1.plot(x, lre, c='r')
    ax1.plot(x, diff, c='b')
    ax1.plot(x, abs_diff, c='g')
leg = ax1.legend()

plt.show()

