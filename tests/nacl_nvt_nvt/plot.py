import numpy as np
import matplotlib.pyplot as plt

headers = {
        "tstep": 0,
        "t": 1,
        "T": 2,
        "P": 3,
        "V": 4,
        "Etot": 5,
        "Eke": 6, 
        "Epe": 7,
        "Epair": 8
        }


data = np.loadtxt("work/trajectory_nvt_5_3.seamm_trj", skiprows=2)
plt.plot(data[:,headers["Etot"]])
plt.show()
