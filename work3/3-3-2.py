import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

delta_t = 0.001

def f(y):
    sigma = 10
    b = 8/3
    r = 28
    return np.array([sigma*(y[1]-y[0]),r*y[0]-y[1]-y[0]*y[2],y[0]*y[1]-b*y[2]])


def next(v):
    k1 = f(v[-1])
    k2 = f(v[-1] + delta_t / 2 * k1)
    k3 = f(v[-1] + delta_t / 2 * k2)
    k4 = f(v[-1] + delta_t * k3)
    next = np.array(v[-1]) + delta_t * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return next


def runge_kutta(v,t):
    for i in range(1, len(t)+1):
        v.append(next(v))


# main
if __name__ == '__main__':
    t = np.arange(0,20,delta_t)

    y_list = []
    q_list = [0., 10 ** (-2), 10 ** (-3)]
    for q in q_list:
        y = []
        y.append(np.array([1., 0., 20.]))
        y[0][2] += q
        runge_kutta(y, t)
        y = np.array(y)
        y_list.append(y)
        # print(y[:, 0])
        # exit()




    # plt.figure()
    # plt.plot(y[0], y[1])

    # fig = plt.figure()
    # ax = Axes3D(fig)
    # plt.plot(y_list[0][::20, 0], y_list[0][::20, 1], y_list[0][::20, 2], c='red')
    # plt.plot(y_list[1][::20, 0], y_list[1][::20, 1], y_list[1][::20, 2], c='blue')
    # plt.plot(y_list[2][::20, 0], y_list[2][::20, 1], y_list[2][::20, 2], c='green')
    # plt.show()

    # gif
    fig = plt.figure()
    ax = Axes3D(fig)
    ims = []
    rate = 200
    for i in range(int(len(t)/rate)):
        line0 = ax.plot(y_list[0][:i * rate:20, 0], y_list[0][:i * rate:20, 1], y_list[0][:i * rate:20, 2], c='red')
        line1 = ax.plot(y_list[1][:i * rate:20, 0], y_list[1][:i * rate:20, 1], y_list[1][:i * rate:20, 2], c='blue')
        line2 = ax.plot(y_list[2][:i * rate:20, 0], y_list[2][:i * rate:20, 1], y_list[2][:i * rate:20, 2], c='green')
        ims.append(line0 + line1 + line2)

    ani = animation.ArtistAnimation(fig, ims)
    ani.save('gif/3-3-2.gif', writer="imagemagick")
    plt.show()