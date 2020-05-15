import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import time


def gauss(traj, s, time): 
    ret = np.sum(traj[:time,3]*np.exp(-(s - traj[:time,1])**2/(2*traj[:time,2]**2)))
    return ret

hills_ref = np.loadtxt('work/output.dat')
start = 10000
stop = 38000
stride = 1000
s_min = 2.4
s_max = 8.0

s = np.linspace(s_min, s_max, 100)

fig = plt.figure()

#creating a subplot
ax1 = fig.add_subplot(1,1,1)

f_s_ref = np.zeros(len(s))

for time in range(start, stop, stride):
    for idx, ss in enumerate(s):
        f_s_ref[idx] -= gauss(hills_ref, ss, time)

f_s_ref /= len(range(start, stop, stride))
f_s_ref -= f_s_ref.min()

ax1.clear()
ax1.plot(s, f_s_ref) 
# this is an inset axes over the main axes
plt.title("Time step " + str(hills_ref[-1, 0]) + "\n width: " + str(hills_ref[-1,2]) + " height: " + str(hills_ref[-1,3]) + "\n cutoff: 5.0 accuracy: 10^-2")
plt.xlabel("Distance (Angstroms)")
plt.ylabel("G (kcal/mol)")

plt.show()
