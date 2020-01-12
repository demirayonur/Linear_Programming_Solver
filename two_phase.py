
"""

This python module implements two phase
approach to given standard Linear Program
in order to find an initial basic feasible
solution when it is not trivial, in other
words, find an initial basic feasible solution
when there is at least one equality or greater
than or equal to constraint.

It first converts the standard Lp into
artificial form then implements PHASE 1.
If PHASE 1 does not result infeasibility, then
return a standard Linear Program to be solved
by the revised simplex algorithm we developed
in the PHASE 2.


Artificial Form:
---------------

min z= 1.x_a
s.t;  Ax + x_a = b
       x,  x_a >= 0

If z* = 0, then continue with PHASE 2.
Otherwise, the problem is infeasible.

"""

import revised_simplex
import numpy as np
import standard_lp
import copy


class TwoPhase(object):

    def __init__(self, lp_model):

        """
        % Construction Function
        :param lp_model: SLP object
        """

        self.lp_model = lp_model
        self.artificial_model = copy.deepcopy(lp_model)  # Artificial model will be based on this attribute
        self.n_x_a = 0  # Number of artificial variables

    def to_artificial(self):

        """
        Finds corresponding artificial form
        of given standard LP

        :return: % INPLACE %
        """

        # First we should find which constraints (rows) requires an artificial variable.
        problematic_rows = []
        for row in range(self.lp_model.m):
            need_artificial_var = True
            for var in range(self.lp_model.n):
                col = self.lp_model.A[:, var]
                n_zeros = (col == 0).sum()
                if n_zeros == self.lp_model.m - 1 and col[row] == 1:
                    need_artificial_var = False
                    break
            if need_artificial_var:
                problematic_rows.append(row)

        # Add artificial columns to A Matrix
        A_prime = self.lp_model.A
        for row in problematic_rows:
            col = np.zeros(self.lp_model.m)
            col[row] = 1
            m, n = A_prime.shape
            actual_A = np.empty(shape=(m, n + 1))
            actual_A[:, :n] = A_prime
            actual_A[:, -1] = col
            A_prime = actual_A
        self.artificial_model.A = A_prime

        # Cost vector of the artificial form
        new_cost_vector = [0 * i for i in range(self.lp_model.n)]
        for i in range(len(problematic_rows)):
            new_cost_vector.append(1)
        new_cost_vector = np.array(new_cost_vector)
        self.artificial_model.c = new_cost_vector

        # Update m & n of the artificial model
        self.artificial_model.m, self.artificial_model.n = self.artificial_model.A.shape

        # Find the number of artificial variables
        self.n_x_a = len(problematic_rows)
        print('Artificial variables are appended to the problem')

    def solve_phase_1(self):

        """
        This function corresponds to phase 1,
        solves the artificial model and if
        feasibility conditions are satisfied,
        then returns the standard LP to be continued
        with the second phase 2 which is expected to
        be solved by revised simplex we developed.

        :return: SLP object
        """

        # Solve artificial model by Revised Simplex Algorithm

        solver = revised_simplex.RevisedSimplex(self.artificial_model)
        sol, obj_val = solver.run()

        # Feasibility Check
        if obj_val != 0.0:
            print('The problem is infeasible')
            return None

        # If there is not any problem with feasibility find phase 2 model
        # Generate new LP to be solved in phase 2 by revised simplex

        x_b = [var for var, val in enumerate(sol) if val > 0]
        x_n = [var for var, val in enumerate(sol) if val == 0 and var < self.lp_model.n]

        A = np.empty(shape=(self.lp_model.m, self.lp_model.n))
        for var in range(self.lp_model.n):
            if var in x_b:
                ind = x_b.index(var)
                col = np.zeros(self.lp_model.m)
                col[ind] = 1
                A[:, var] = col
            else:
                A[:, var] = solver.B_inv[:, x_n.index(var)]
        phase2_model = standard_lp.SLP(self.lp_model.c, A, [sol[var] for var in x_b], mode=self.lp_model.mode)
        print('Phase 1 is completed')
        return phase2_model
