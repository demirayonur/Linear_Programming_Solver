"""

This module includes RevisedSimplex
class which is responsible for finding
optimal solution given standard LP model

"""

import numpy as np
import utils_


class RevisedSimplex(object):

    def __init__(self, lp_model):

        """
        % Construction Function

        :param lp_model: SLP object
        """

        self.lp_model = lp_model
        self.x_b, self.x_n = [], []
        self.c_b, self.c_n = None, None
        self.B_inv = None
        self.N = None

        self.get_trivial_initial_bfs()

    def get_trivial_initial_bfs(self):

        """
        This function finds an initial basic
        feasible solution to the system when
        it is trivial, the case when all
        constraints are less than or equal to
        type constraints, and updates 'x_b',
        'x_n', 'B_inv', 'N', 'c_b', 'c_n'

        :return: % INPLACE %
        """

        # First determine basic variables
        for row in range(self.lp_model.m):
            # for given row, we should find which variable or column of A satisfies canonical form
            for var in range(self.lp_model.n):
                # it should include 1 in the 'row' and 0 for the others
                col = self.lp_model.A[:, var]
                num_zeros = (col == 0).sum()
                if num_zeros == self.lp_model.m - 1 and col[row] == 1:
                    self.x_b.append(var)
                    break

        # Then non-basic variables (if a variable does not take place in 'x_n')
        self.x_n = [i for i in range(self.lp_model.n) if i not in self.x_b]

        # Initial 'B_inv' is an identity matrix with m x m
        self.B_inv = np.identity(self.lp_model.m)

        # Cost of variables in 'x_b' constitutes 'c_b'
        self.c_b = np.array([self.lp_model.c[i] for i in self.x_b])

        # Cost of variables in 'x_n' constitutes 'c_n'
        self.c_n = np.array([self.lp_model.c[i] for i in self.x_n])

        # N
        self.N = np.empty((self.lp_model.m, self.lp_model.n - self.lp_model.m))
        for j in range(self.N.shape[1]):
            self.N[:, j] = self.lp_model.A[:, self.x_n[j]]

    def run(self):

        """
        Iterates revised simplex algorithm until
        it finds the optimal solution and returns
        the value of x and objective function value
        respectively.
        :return: x, obj_val
        """
        x = []
        obj_val = 0
        iter_no = 0
        while True:  # It will be broken when the optimality condition or unboundedness conditions are met.

            print('iteration no: ', str(iter_no))
            # Compute Parameters
            w = np.matmul(self.c_b, self.B_inv)
            b_bar = np.matmul(self.B_inv, self.lp_model.b)
            reduced_costs = np.matmul(w, self.N) - self.c_n

            print('\tobjective function value: ', str(np.matmul(w, self.lp_model.b)))

            # Optimality Check
            is_optimal = utils_.optimality_check(reduced_costs)
            if is_optimal:
                for var in range(self.lp_model.n):
                    if var in self.x_b:
                        ind = self.x_b.index(var)
                        x.append(b_bar[ind])
                    else:
                        x.append(0)
                obj_val = np.matmul(w, self.lp_model.b)
                break

            # Entering Variable
            k_ind = int(np.argmax(reduced_costs))  # index of entering variable
            k_var = self.x_n[k_ind]  # variable name of entering variable
            k_cost = self.c_n[k_ind]  # cost of this entering variable
            print('\tentering variable: ', str(k_var))

            y_k = np.matmul(self.B_inv, self.lp_model.A[:, k_var])

            # Unboundedness Check
            is_unbounded = utils_.unboundedness_check(y_k)
            if is_unbounded:
                break

            # Leaving Variable
            r_ind = utils_.min_ratio_test(self.lp_model.m, y_k, b_bar)  # index of leaving variable
            r_var = self.x_b[r_ind]  # variable name of leaving variable
            r_cost = self.c_b[r_ind]  # cost of this leaving variable
            print('\tleaving variable: ', str(r_var))

            # Update Solution
            self.x_b[r_ind] = k_var
            self.x_n[k_ind] = r_var
            self.c_b[r_ind] = k_cost
            self.c_n[k_ind] = r_cost
            self.N[:, k_ind] = self.lp_model.A[:, r_var]

            # Update B inverse
            E = np.identity(self.lp_model.m)
            col = np.empty(shape=self.lp_model.m)
            for i in range(self.lp_model.m):
                if i == r_ind:
                    col[i] = 1.0 / y_k[i]
                else:
                    col[i] = -y_k[i] / y_k[r_ind]
            E[:, r_ind] = col
            self.B_inv = np.matmul(E, self.B_inv)

            # Update Iteration Number
            iter_no += 1

        return x, obj_val
