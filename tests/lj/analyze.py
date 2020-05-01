import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import time

kb = 0.0019872041
T = 130

os.system('./tcp.sh &')
time.sleep(3)

fig = plt.figure()

#creating a subplot
ax1 = fig.add_subplot(1,1,1)

def animate(i):
    hills = np.loadtxt('work/output.dat')

    start = 0
    stop = len(hills) 
    stride = 10
    x_min = 1.0
    x_max = 14.0

    ax1.clear()
    x = np.linspace(x_min, x_max, 100)
    f_x = np.zeros(len(x))
    
    gauss = lambda x, n: np.sum(hills[:n,3]*np.exp(-(x - hills[:n,1])**2/(2*hills[:n,2]**2)))

    for n in range(start, stop, stride):
        f_x -= np.array([gauss(xx, n) for xx in x])

    f_x /= len(range(start, stop, stride))
    f_x -= f_x.min()

    epsilon = 0.2
    sigma = 3.88
    f_ref = 4*epsilon * ((sigma/x)**12 - (sigma/x)**6) + 2*kb*T*np.log(x)
    f_ref -= f_ref.min()

    ax1.clear()
    ax1.plot(x, f_x, x, f_ref)

    ax1.set_ylim([0.0, 1.75])
    ax1.set_xlim([3.5, 15])

    ax1.set_label('Distance')
    #plt.label('Free energy')
    plt.title('Time step: '+ str(hills[-1, 0]))

    # this is an inset axes over the main axes

    axin1 = ax1.inset_axes([0.1, 0.6, 0.3, 0.3])
    axin1.plot(hills[-100:,0], hills[-100:,1])

ani = animation.FuncAnimation(fig, animate, interval=10)
plt.show()
