import numpy as np
import matplotlib.pyplot as plt

with open ("input.dat") as f:
    a = f.readlines()

nx = int(a[0].strip().split()[2])
nz = int(a[1].strip().split()[2])
cnftime = float(a[5].strip().split()[2])
simtime = float(a[6].strip().split()[2])

times = np.arange(cnftime, simtime+1e-10, cnftime)

N = nx * nx * nz
print(N)

x0, y0, z0 = [], [], []
msd = np.array([0])

with open("cnf.000") as f:
    a = f.readlines()
    for i in range(2, N+2):
        x0.append(float(a[i].strip().split()[0]))
        y0.append(float(a[i].strip().split()[1]))
        z0.append(float(a[i].strip().split()[2]))


for cnftag in range(1, len(times)+1):
    x1, y1, z1 = [], [], []
    temp = 0
    with open("cnf." + "{0:0=3d}".format(cnftag)) as f:
        a = f.readlines()
        for i, j  in zip(range(2, N+2), range(N)):
            x1.append(float(a[i].strip().split()[0]))
            y1.append(float(a[i].strip().split()[1]))
            z1.append(float(a[i].strip().split()[2]))
            temp += pow(x1[j] - x0[j], 2)
        print(temp/N)
        msd = np.append(msd, float(temp/N))

times = np.insert(times, 0, 0)
plt.plot(times, msd)
plt.show()

