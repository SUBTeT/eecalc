import numpy as np

# main
if __name__ == '__main__':
    s = 0.03
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

    S_t = np.linalg.matrix_power(S, 60)
    g_next = np.dot(S_t, g)

    print(S_t)