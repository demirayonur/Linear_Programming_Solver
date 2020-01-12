import standard_lp
import numpy as np


# Example 7.2.1 in Metin Türkay's Lecture Notes

c = np.array([-19, -13, -12, -17, 0, 0, 0])  # cost vector
b = np.array([225, 117, 420])  # RHS vector
A = np.array([[3, 2, 1, 2, 1, 0, 0], [1, 1, 1, 1, 0, 1, 0], [4, 3, 3, 4, 0, 0, 1]])  # A matrix
mode = 'min'  # in order to maximize the objective function, you should assign 'max' to 'mode' variable

model = standard_lp.SLP(c, A, b, mode)  # create LP object
model.solve()  # solve LP by Revised Simplex

print('\noptimal solution is: ', str(model.get_variable_values()))  # Get optimal solution
print('\noptimal objective function value is: ', str(model.get_obj_value()))  # Get optimal objective function value
