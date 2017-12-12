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

with open("Energy_graph.txt") as f:
    lines = f.readlines()
    for line in lines:
        if line==lines[0]:
            if x:
                ax1.plot(x, lre, c='r', label='energy curve')
                ax1.plot(x, diff, c='b', label='difference energy curve')
            x = []
            lre = []
            diff = []
        else : 
            try:
                l = line.split(";")
       
                x.append(l[0])
                lre.append(l[1])
                diff.append(l[2])
            except:
               print(line)




leg = ax1.legend()

plt.show()

#a la fin t'es obligé de le supprimer ou de tout effacer sinon java réecrit a la suite
#os.remove("Energy_graph.txt")
