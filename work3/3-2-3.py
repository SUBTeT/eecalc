import numpy as np
from numba import jit
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def f(t, x):
    return np.array([x[1], -x[0]])


def next(t, v):
    k1 = f(None, v[-1])
    k2 = f(None, v[-1] + delta_t * k1)
    next = np.array(v[-1]) + delta_t * (k1 + k2) / 2
    return next


def heun(v, stop):
    for i in range(1, int(stop)):
        v.append(next(None, v))


# main
if __name__ == '__main__':
    tau = 2 * np.pi
    delta_t = tau / 64

    v_c = []
    v_c.append(np.array([0, 1]))
    heun(v_c, 5 * tau / delta_t)

    t = []
    for i in range(int(5 * tau / delta_t)):
        t.append(delta_t * i)

    v_a = []
    for i in range(int(5 * tau / delta_t)):
        v_a.append([np.sin(t[i]), np.cos(t[i])])

    error = []
    for i in range(int(5 * tau / delta_t)):
        error.append(np.linalg.norm(v_c[i] - v_a[i]))

    plt.figure()
    plt.plot(t, error)

    ##############
    delta_t_list = []
    p_list = range(3, 19)
    for p in p_list:
        delta_t_list.append(tau * 2 ** (-p))

    log2_E_r = []

    for delta_t in delta_t_list:
        v_c = []
        v_c.append(np.array([0, 1]))
        heun(v_c, 5 * tau / delta_t)

        t = []
        for i in range(int(tau / delta_t)):
            t.append(delta_t * i)

        v_a = []
        for i in range(int(tau / delta_t)):
            v_a.append([np.sin(t[i]), np.cos(t[i])])

        error = []
        for i in range(int(tau / delta_t)):
            error.append(np.linalg.norm(v_c[i] - v_a[i]))

        log2_E_r.append(np.log2(max(error)))

    plt.figure()
    plt.plot(p_list, log2_E_r)
    a, b = np.polyfit(p_list, log2_E_r, 1)
    print(a)
    #############
    plt.figure()
    plt.scatter(np.array(v_c).T[0], np.array(v_c).T[1])
    plt.show()


    '''
    # gif
    fig = plt.figure()

    ims = []
    for i in range(int(5 * tau / delta_t / 3)):
        dot = plt.scatter(v_c[i * 3][0], v_c[i * 3][1], c='red')
        ims.append([dot])

    ani = animation.ArtistAnimation(fig, ims)
    ani.save('gif/3-2-3.gif', writer="imagemagick")
    plt.show()
    '''
