import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import time


def gauss(traj, s, time): 
    ret = np.sum(traj[:time,3]*np.exp(-(s - traj[:time,1])**2/(2*traj[:time,2]**2)))
    return ret

hills_ssages = np.loadtxt('data/hills_ref.dat')

hills_ssages[:,1] *= 10.0
hills_ssages[:,2] *= 10.0
#hills_ssages[:,3] /= 4.184

# SSAGES

start = 10000
stop = 50000
stride = 1000
s_min = 1.0
s_max = 7.0

s_ssages = np.linspace(s_min, s_max, 100)

fig = plt.figure()

#creating a subplot
ax1 = fig.add_subplot(1,1,1)

f_s_ssages = np.zeros(len(s_ssages))

for time in range(start, stop, stride):
    for idx, ss in enumerate(s_ssages):
        f_s_ssages[idx] -= gauss(hills_ssages, ss, time)

f_s_ssages/= len(range(start, stop, stride))
f_s_ssages-= f_s_ssages.min()

ax1.clear()
#ax1.plot(s_ssages, f_s_ssages) 

# MDI

hills_ref = np.loadtxt('work/output.dat')

hills_ref[:,3] *= 4.184

start_mdi = 10000
stop_mdi = 50000
stride_mdi = 1000
s_min_mdi = 1.0
s_max_mdi = 7.0

s_mdi = np.linspace(s_min_mdi, s_max_mdi, 100)

f_s_mdi = np.zeros(len(s_mdi))

for time in range(start_mdi, stop_mdi, stride_mdi):
    for idx, ss in enumerate(s_mdi):
        f_s_mdi[idx] -= gauss(hills_ref, ss, time)

f_s_mdi /= len(range(start_mdi, stop_mdi, stride_mdi))
f_s_mdi -= f_s_mdi.min()

# Paper

data = np.loadtxt("data/g_ssages.dat")
data[:,0] *= 10

ax1.plot(s_mdi, f_s_mdi, color='r') 
ax1.plot(data[:,0], data[:, 1], color='m') 
# this is an inset axes over the main axes
plt.title("Time step " + str(hills_ref[-1, 0]) + "\n width: " + str(hills_ref[-1,2]) + " height: " + str(hills_ref[-1,3]) + "\n cutoff: 10.0 accuracy: 10^-5")
ax1.set_ylim([0,20])
ax1.set_xlim([1,7])
plt.xlabel("Distance (Angstroms)")
plt.ylabel("G (kcal/mol)")

plt.show()
