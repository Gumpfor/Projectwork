import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys


filename = str(sys.argv[1]).zfill(3) + ".x"
with open(filename) as f:
    a = f.readline().strip().split()
    rc = [float(a[i]) for i in range(3)]
print(rc)


filename = str(sys.argv[1]).zfill(3) + ".plt"
with open(filename) as f:
    temp = f.readline()
    for dat in f.readlines():
        line = dat.strip().split()
        r = [float(line[i]) for i in range(3)]
        r = np.matrix(r)
        mat = np.asarray([float(line[i]) for i in range(3, 12)])
        mat = np.matrix(mat.reshape((3,-1)))
        print(np.dot(np.dot(rc-r, mat), (rc-r).transpose()))

print("\n")
