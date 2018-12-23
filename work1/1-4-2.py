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

    w, v = np.linalg.eig(S)
    print(w)
    print(v)