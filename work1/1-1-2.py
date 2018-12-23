import matplotlib.pyplot as plt
import numpy as np

# equation
def func_f(x):
    return (5*x**3-3*x)/2

# Newton's method
def newton(func_f, x, i):
    num_calc = 0 # the number of calculations
    x_tru = np.sqrt(3 / 5)  # true solution
    y = []
    for j in range(0, i):
        x_next = x - func_f(x) / ((15 * x ** 2 - 3) / 2)
        dif = abs(x_next - x_tru)

        num_calc += 1
        x = x_next
        y.append(np.log2(dif))

    print("x = {}".format(x))
    plt.xlabel("number of trials")
    plt.ylabel("log_2(error)")
    plt.plot(y)
    plt.show()
    return x

# main
if __name__ == '__main__':
    # calculate solution (equation, the first approximate solution, the number of calculations）
    solution = newton(func_f, 0.5, 9)