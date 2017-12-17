import numpy as np
import os
import matplotlib.pyplot as plt

first = True
colours = ['r', 'b', 'g', 'gray', 'pink', 'orange']
def plot_line(x, y, labels):
    global first
    if first:
        first = False
        for (i, yi) in enumerate(y):
            ax1.plot(x[i], yi, c = colours[i%len(colours)], label=labels[i])
    else:
        for (i, yi) in enumerate(y):
            ax1.plot(x[i], yi,  c = colours[i%len(colours)])


fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Running time")
ax1.set_xlabel('Iterations')
ax1.set_ylabel('Running time')
x = []
xx = []
time = []
times = []
n = []
names = []
with open("running_time.txt") as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith("NEW"):
            l = line.split('/')
            if x:
                xx.append(x)
                times.append(time)
            names.append(l[-1][:-5])
            x, time = [], []
        else : 
            try:
                l = line.split(";")
                if len(l)==2:
                    x.append(float(l[0]))
                    time.append(float(l[1]))
                else:
                    n.append(l[0])
                    print(l)
            except:
               print(line)
               
              
xx.append(x)
times.append(time)
plot_line(xx, times, names)
leg = ax1.legend()
#plt.show()

times_trunc = [time[50:] for time in times]

plt.figure()
n2 = [int(num) for num in n]
amp = max(n2)-min(n2)
plt.boxplot(times_trunc, 1, '', positions=n2, widths = amp//20)
plt.xlim(min(n2)-amp//20, max(n2)+amp//20)
plt.show()

