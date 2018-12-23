import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 100
delta_x = 1.0 / N
delta_t = 0.001
c = 5
T_D = 0.1
mu = c ** 2 * delta_t ** 2 / delta_x ** 2

limit = T_D / delta_t * 4  # period

x = []
for i in range(N):
    x.append(i * delta_x)

t = []
for i in range(int(limit)):  # check
    t.append(i * delta_t)


def calc_next(u):
    u_k_next = np.zeros(N)

    # fixed end
    u_k_next[N - 1] = -u[-2][N - 1] + mu * u[-1][N - 2] + 2 * (1 - mu) * u[-1][N - 1]
    # free end
    # u_k_next[N - 1] = -u[-2][N - 1] + mu * u[-1][N - 2] + (2 - mu) * u[-1][N - 1]

    if len(u) * delta_t < T_D:
        u_k_next[0] = 1.

    for i in range(1, N - 1):
        u_k_next[i] = -u[-2][i] + mu * u[-1][i - 1] + 2 * (1 - mu) * u[-1][i] + mu * u[-1][i + 1]

    return u_k_next


def make_u(u):
    for i in range(int(limit)):
        u_k_next = calc_next(u)
        u.append(u_k_next)
    return u


# main
if __name__ == '__main__':
    # initial
    u = []
    a = np.zeros(N)
    a[0] = 1
    u.append(a)
    u.append(a)

    make_u(u)

    # print(u[2500])
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.plot(u[100])
    # plt.show()
    #
    # print(len(u))
    # exit()


    # gif
    fig = plt.figure()

    ims = []
    for i in range(int(len(u)/10)):
        line, = plt.plot(u[i * 10], "r")
        ims.append([line])

    ani = animation.ArtistAnimation(fig, ims)
    ani.save('gif/3-1-2.gif', writer="imagemagick")
    plt.show()
