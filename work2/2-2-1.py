import numpy as np
from numba import jit

error = 1e-10
limit = 1e+7


@jit
def calc_next(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> np.ndarray:
    n = len(b)
    x_next = np.zeros(n)
    for i in range(0, n):
        sum = 0
        for j in range(0, i):
            sum += A[i][j] * x[j]
        for j in range(i + 1, n):
            sum += A[i][j] * x[j]
        x_next[i] = (b[i] - sum) / A[i][i]
    return x_next


def solve_by_jacobi(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> np.ndarray:
    while True:
        x_next = calc_next(A, b, x)
        dif = np.abs(np.linalg.norm(x_next) - np.linalg.norm(x))
        if dif <= error:
            break
        if dif > limit:
            print("Jacobi method cannot get the solution...")
            exit(1)
        x = x_next
    return x_next


if __name__ == '__main__':
    A = np.array([[3.0, 2.0, 1.0], [1.0, 3.0, -2.0], [2.0, -1.0, 4.0]])
    b = np.array([4, 6, -3])

    # initial solution
    x = np.zeros(len(b))

    print("A: \n", A)
    print("b: \n", b)
    print("x: \n", solve_by_jacobi(A, b, x))
