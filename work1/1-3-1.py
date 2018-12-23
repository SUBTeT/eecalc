import numpy as np
import matplotlib.pyplot as plt

# main
if __name__ == '__main__':
    s = 0.03 # initial settings
    p = q = 0.485
    g = np.array([1,0,0,0,0,0,0])
    S = np.array([[s+p,p,0,0,0,0,0],
                  [q,s,p,0,0,0,0],
                  [0,q,s,p,0,0,0],
                  [0,0,q,s,p,0,0],
                  [0,0,0,q,s,p,0],
                  [0,0,0,0,q,s,p],
                  [0,0,0,0,0,q,s+q]
                  ])

    y = [] # list of the results
    for i in range(0,7):
        y.append([])

    for t in range(0, 60):
        for i in range(0,7): # record the 7 results
            y[i].append(g[i])
        g_next = np.dot(S, g)
        g = g_next # refresh g

    print(g)
    ptn = ["b", "g", "r", "c", "m", "y", "k"]
    for i in range(0,7):
        plt.plot(y[i], ptn[i],label="Point{}".format(i))
    plt.legend()
    plt.xlabel("time[s]")
    plt.ylabel("probability")
    plt.show()