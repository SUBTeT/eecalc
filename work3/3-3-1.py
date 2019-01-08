import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

delta_t = 0.001

def f(t, y):
    sigma = 10
    b = 8/3
    r = 28
    return np.array([sigma*(y[1]-y[0]),r*y[0]-y[1]-y[0]*y[2],y[0]*y[1]-b*y[2]])


def next(t, v):
    k1 = f(None, v[-1])
    k2 = f(None, v[-1] + delta_t / 2 * k1)
    k3 = f(None, v[-1] + delta_t / 2 * k2)
    k4 = f(None, v[-1] + delta_t * k3)
    next = np.array(v[-1]) + delta_t * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return next


def runge_kutta(v,t):
    for i in range(1, len(t)+1):
        v.append(next(None, v))


# main
if __name__ == '__main__':
    t = np.arange(0,10,delta_t)

    y = []
    y.append(np.array([1, 0, 20]))
    # print(v[0][1])
    runge_kutta(y,t)
    y = np.array(y)
    # print(y[:, 0])
    # exit()


    # plt.figure()
    # plt.plot(y[0], y[1])

    # fig = plt.figure()
    # ax = Axes3D(fig)
    # dot = plt.plot(y[::20, 0], y[::20, 1], y[::20, 2], c='red')
    # plt.show()

    # gif
    fig = plt.figure()
    ax = Axes3D(fig)


    ims = []
    for i in range(int(len(t)/200)):
        line = plt.plot(y[:i*200:20, 0], y[:i*200:20, 1],y[:i*200:20, 2], c='red')
        ims.extend([line])

    ani = animation.ArtistAnimation(fig, ims)
    ani.save('gif/3-3-1.gif', writer="imagemagick")
    plt.show()
