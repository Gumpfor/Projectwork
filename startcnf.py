from math import pi
import sys
import copy
from random import *

with open("input.dat") as f:
    a = f.readlines()

nx  = int(a[0].strip().split()[2])
nz  = int(a[1].strip().split()[2])
box = float(a[2].strip().split()[2])
ha  = float(a[3].strip().split()[2])
hb  = float(a[4].strip().split()[2])
cnftime  = a[5].strip().split()[2]
simtime  = a[6].strip().split()[2]

null = 0.
n = nx * nx * nz

distz = (box - 2*hb*nz) / nz
distx = (box - 2*ha*nx) / nx

if distx < 0 or distz < 0:
    sys.exit("Overlap!!!!")

koord = [-box/2+distx/2+ha, -box/2+distx/2+ha, -box/2+distz/2+hb]

# print("nx = ", nx, "\nnz = ", nz, "\nbox = ", box, "\nha = ", ha, "\nhb = ", hb)

f = open("inp.cnf", "w")

f.write(str(n) + "\n")
f.write(str(box) + "\n")
f.write("0.0\n")
f.write("rx          ry          rz          zax        zay        zaz        xax        xay        xaz        vx          vy          vz          omegax      omegay      omegaz\n")

for i in range(nx):
    for j in range(nx):
        for k in range(nz):
            a = copy.deepcopy(koord)
            a[0] += i*(distx+2*ha)
            a[1] += j*(distx+2*ha)
            a[2] += k*(distz+2*hb)
            f.write(("%+9.7f  %+9.7f  %+9.7f  %9.7f  %9.7f  %9.7f  %9.7f  %9.7f  %9.7f  %+9.7f  %+9.7f  "
                    "%+9.7f  %+9.7f  %+9.7f  %+9.7f\n") %(a[0], a[1], a[2],
                    null, null, hb, ha, null, null,
                    uniform(-1, 1), uniform(-1, 1), uniform(-1, 1),
                    uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)))
            

print(("Raumfuellungsfaktor:             %6.5f") %(n*ha*hb*pi/pow(box, 3)))

fi = open("input", "w")
fi.write(str(ha) + "\n" + str(hb) + "\n" + cnftime + "\n" + simtime + "\n")
fi.close()
f.close()
