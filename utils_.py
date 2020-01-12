"""

This python module does not include any class.
However, it consists of functions which are
necessary to implement the algorithms.

"""

import numpy as np
import math


def optimality_check(reduced_costs):

    """
    Given reduced cost array, checks whether
    optimality condition are satisfied or not.
    If satisfied returns True; otherwise False
    :param reduced_costs: 1-D Numpy Array
    :return: Boolean
    """

    flag = False
    if np.max(reduced_costs) <= 0:
        print('Optimal solution is found')
        flag = True

    return flag


def unboundedness_check(y_k):

    """
    Given y_k vector, checks whether unboundedness
    conditions are satisfied or not. If satisfied
    returns True; otherwise returns False
    :param y_k: 1-D Numpy Array
    :return: Boolean
    """

    flag = False
    if np.max(y_k) <= 0:
        print('The model is unbounded')
        flag = True
    return flag


def min_ratio_test(m, y_k, b_bar):

    """
    :param m: Integer
    :param y_k: 1-D Numpy Array
    :param b_bar: 1-D Numpy Array
    :return: int, index of the leaving variable in 'x_b'
    """

    r, temp = -1, math.inf
    for i in range(m):
        if y_k[i] > 0 and b_bar[i] / y_k[i] < temp:
            r, temp = i, b_bar[i] / y_k[i]

    return r
