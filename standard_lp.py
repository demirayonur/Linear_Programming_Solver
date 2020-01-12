"""
This python module consists of
SLP class which is abbreviation
of 'Standard Linear Program'.This
class defines the following parts
of any Standard Linear Program:

--> c: Cost vector
--> b: Right hand side vector
--> A: Technology coefficient matrix
--> mode: 'min' or 'max'

Caution: The algorithms are designed
to minimize the objective function.
Therefore, if the mode is 'max' then
cost vector will be multiplied with
minus 1 and the rest will be the same.
"""

import revised_simplex
import numpy as np
import two_phase


class SLP(object):

    def __init__(self, c, A, b, mode):

        """
        The standard form is as follows:
        min cx
        st; Ax = b

        :param c: (1,n)
        :param A: (m,n)
        :param b: (m,1)
        :param mode: 'min' for minimization
                     'max' for maximization
        """

        # Parameters
        self.m, self.n = A.shape  # m: number of constraints, n: number of decision variables
        self.c = c
        self.A = A
        self.b = b

        # Maximization or Minimization
        self.mode = mode
        if self.mode == 'max':
            self.c = -1 * c

        # Decision Variables and Objective Function Value (will be filled after solved by Revised Simplex)
        self.x = []
        self.objective_value = None

    def get_variable_values(self):

        """
        Returns the decision variable values
        as a 1-D Numpy Array

        :return: A 1-D Numpy Array
        """

        return np.array(self.x)

    def get_obj_value(self):

        """
        Returns the optimal objective
        function value

        :return: Float
        """

        return self.objective_value

    def solve(self):

        # First determine we apply 2-Phase method or not before continue with revised simplex.

        problematic_rows = []
        for row in range(self.m):
            need_artificial_var = True
            for var in range(self.n):
                col = self.A[:, var]
                n_zeros = (col == 0).sum()
                if n_zeros == self.m - 1 and col[row] == 1:
                    need_artificial_var = False
                    break
            if need_artificial_var:
                problematic_rows.append(row)

        # When the initial BFS is trivial
        if len(problematic_rows) == 0:
            solver = revised_simplex.RevisedSimplex(self)

        # When the initial BFS is not trivial
        else:
            two_phase_solver = two_phase.TwoPhase(self)
            two_phase_solver.to_artificial()
            phase2_model = two_phase_solver.solve_phase_1()
            solver = revised_simplex.RevisedSimplex(phase2_model)

        # Solve
        self.x, self.objective_value = solver.run()
