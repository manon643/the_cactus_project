import numpy as np
import os
import matplotlib.pyplot as plt

with open("Energy_data.txt") as f:
    lines = f.readlines()
    x = [line.split()[0] for line in lines]
    y = [line.split()[1] for line in lines]


fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("Energy graph")
ax1.set_xlabel('Iterations')
ax1.set_ylabel('Energy')

ax1.plot(x,y, c='r', label='energy curve')

leg = ax1.legend()

plt.show()

#a la fin t'es obligé de le supprimer ou de tout effacer sinon java réecrit a la suite
os.remove("Energy_graph.txt")