import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys


def coordinates(A, center):
    # find the rotation matrix and radii of the axes
    U, s, rotation = linalg.svd(A)
    radii = 1.0/np.sqrt(s)

    # now carry on with EOL's answer
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(x)):
        for j in range(len(x)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center
    return x, y, z

def calcmatrix(ax):
    # xaxis /= np.linalg.norm(xaxis)
    # yaxis /= np.linalg.norm(yaxis)
    # zaxis /= np.linalg.norm(zaxis)
    ev = [1,1,0.25]
    V = np.matrix(ax).transpose()
    Vinv = linalg.inv(V)
    L = np.matrix([[ev[0], 0, 0], [0, ev[1], 0], [0, 0, ev[2]]])
    filename = str(sys.argv[1]).zfill(3) + ".mat"
    # print(
    # f = open(filename, "a")

    # f.write(str(np.dot(np.dot(V, L), Vinv)))
    return np.dot(np.dot(V, L), Vinv)


def main():
    r = []
    matrices = []
    ell1, ell2 = 0, 0

    filename = str(sys.argv[1]).zfill(3) + ".plt"
    print(filename)
    with open(filename) as f:
        temp = f.readline()
        try:
            ell1, ell2 = int(temp.strip().split()[0]), int(temp.strip().split()[1])
        except:
            pass
        for line in f:
            temp = line.strip().split()
            r.append([float(temp[0]), float(temp[1]), float(temp[2])])
            xax = [float(temp[3]), float(temp[4]), float(temp[5])]
            yax = [float(temp[6]), float(temp[7]), float(temp[8])]
            zax = [float(temp[9]), float(temp[10]), float(temp[11])]
            matrices.append([xax, yax, zax])


    # your ellispsoid and center in matrix form
    # A = np.array([[1,0,0],[0,1,0],[0,0,0.25]])
    # center = [1,0,0]
    # B = np.array([[1,0,0],[0,1,0],[0,0,0.25]])
    # center2 = [3,0,0]

    X, Y, Z = [], [], []

    # print(np.reshape(np.asarray(matrices[0]),[3,3]))

    for i in range(len(r)):
        x, y, z = coordinates(np.reshape(np.asarray(matrices[i]),[3,3]), r[i])
        X.append(x)
        Y.append(y)
        Z.append(z)

    # x, y, z = coordinates(A, center)
    # x1, y1, z1 = coordinates(B, center2)

    # X = [x, x1]
    # Y = [y, y1]
    # Z = [z, z1]

    # plot
    fig = plt.figure(figsize=(10,10))

    def makeplot(position, elev, angle):
        ax = fig.add_subplot(position, projection='3d')
        count = 1
        for px, py, pz in zip(X, Y, Z):
            if count == ell1 or count == ell2:
                ax.plot_wireframe(px, py, pz, rstride=4, cstride=4, color='r', alpha=0.2)
                filename = str(sys.argv[1]).zfill(3) + ".x"
                with open(filename) as f:
                    temprc = [float(tempel) for tempel in f.readline().strip().split()]
                    ax.plot([temprc[0]], [temprc[1]], [temprc[2]], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=5, alpha=0.6)
            else:
                ax.plot_wireframe(px, py, pz, rstride=4, cstride=4, color='b', alpha=0.2)
            count += 1
        # ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.2)
        # ax.plot_wireframe(x1, y1, z1, rstride=4, cstride=4, color='r', alpha=0.2)
        plt.xlabel("x")
        plt.ylabel("y")
        ax.set_xlim3d(-2.5,2.5)
        ax.set_ylim3d(-2.5,2.5)
        ax.set_zlim3d(-2.5,2.5)
        ax.view_init(elev, angle)
        return ax

    ax1 = makeplot(221, 0, 0)
    ax2 = makeplot(222, 0, 90)
    ax2 = makeplot(223, 90, 0)
    ax2 = makeplot(224, 60, 30)
        # Create cubic bounding box to simulate equal aspect ratio
    # max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max()
    # Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(x.max()+x.min())
    # Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(y.max()+y.min())
    # Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(z.max()+z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    # for xb, yb, zb in zip(Xb, Yb, Zb):
    #    ax.plot([xb], [yb], [zb], 'w')
    # plt.show()
    # plt.close(fig)
    # plt.show()
    plt.savefig(str(sys.argv[1]).zfill(3) + ".png")
    del fig


if __name__== "__main__":
    main()
