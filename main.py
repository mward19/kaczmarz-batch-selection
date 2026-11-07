from typing import List
import os
from collections import defaultdict

import numpy as np

from matplotlib import pyplot as plt

from tqdm.auto import tqdm

def get_system(rng, n: int, d: int):
    """Create A x = b system. Returns A, x, and b."""
    A = rng.normal(size=(n, d))
    A /= np.linalg.norm(A, axis=1, keepdims=True)  # Normalize rows
    x = rng.normal(size=d)
    b = A @ x
    return A, x, b

def get_worst_system(s=10):
    A = np.array([
        [1, 0],
        [1, 0],
        [0, 1],
        [0, 1],
        [1, 0],
    ]).astype(float)
    x = np.zeros(2)
    b = np.zeros(5)
    b[-1] = 2 * s

    teacher = np.full(2, s)

    return A, x, b, teacher


class KaczmarzSolver:
    def __init__(self, x0, A, b, x_sol, *, maxiter=100, batch_size=None, seed=None):
        self.x0 = x0
        self.x = self.x0.copy()
        self.k = 0 # Current step index
        self.A = A
        self.n, self.d = A.shape
        self.b = b
        self.x_sol = x_sol
        self.e = self.A @ self.x_sol - self.b # error
        self.max_e = np.max(np.abs(self.e))
        self.maxiter = maxiter
        self.batch_size = self.n if batch_size is None else batch_size
        self.batch_indices = None

        self.true_residuals = self.A @ self.x_sol - self.b

        self.rng = np.random.default_rng(seed=seed)

        self.scores = []
        self.parameter_errors = []
        self.diagnostics = defaultdict(list)

        self.name = type(self).__name__

    def reset(self):
        self.x = self.x0.copy()
        self.k = 0
        self.scores = []
        self.parameter_errors = []
        self.diagnostics = defaultdict(list)
        self._score()

    def get_batch_indices(self):
        # Uniformly sampled
        self.batch_indices = self.rng.choice(np.arange(self.n), size=self.batch_size, replace=False)
        return self.batch_indices

    def _select_index(self):
        indices = self.get_batch_indices()
        return self.select_index(indices)
    
    def select_index(self, indices):
        """Selects the index of the row of A to step on at a given iteration"""
        raise NotImplementedError()
    
    def _step(self):
        i = self._select_index()
        A_row = self.A[i]
        self.x += (self.b[i] - np.inner(A_row, self.x)) / np.sum(A_row**2) * A_row
        self.k += 1

    def converged(self):
        return self.k >= self.maxiter
    
    def _score(self):
        score = np.sum((self.x - self.x_sol)**2)
        self.scores.append(score)
        self.parameter_errors.append(score)

    def calc_diagnostics(self):
        return None
    
    def solve(self):
        if not self.scores:
            self._score()
        while not self.converged():
            print(f'Iterating (k = {self.k})', end='\r')
            self.calc_diagnostics()
            self._step()
            self._score()
        
        return self.x

    def get_subset_indices(self, subset="all"):
        if subset == "all":
            return np.arange(self.n)
        else:
            raise NotImplementedError(f"Index subset {subset} not implemented for {self.name}")


class RandomizedKaczmarz(KaczmarzSolver):
    # uniformly random selection
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.min_sing_val_sq = np.linalg.svdvals(self.A)[-1]**2

    def select_index(self, indices):
        return self.rng.choice(indices)
    
    def calc_diagnostics(self):
        # Expected convergence rate for randomized Kaczmarz (rows are normalized,
        # so the squared Frobenius norm is just the row count)
        bound = (1 - self.min_sing_val_sq / self.A.shape[0]) * np.sum((self.x - self.x_sol)**2)
        self.diagnostics["bounds"].append(bound)
    


class Motzkins(KaczmarzSolver):
    def select_index(self, indices):
        residuals = self.A[indices] @ self.x - self.b[indices]
        return indices[np.argmax(residuals**2)]

    def calc_diagnostics(self):
        # For the upcoming iteration. Lemma 1
        max_residual = np.max(np.abs(self.A @ self.x - self.b))
        if max_residual > 4 * self.max_e: # not close to sol'n
            bound = np.sum((self.x - self.x_sol)**2) - 1/2 * max_residual**2
        else:
            bound = float('nan') # temporary placeholder
        self.diagnostics["bounds"].append(bound)


class TeacherOpt(KaczmarzSolver):
    def __init__(self, *args, teacher, noised_indices=[], bound_type=0, **kwargs):
        self.teacher = teacher

        super().__init__(*args, **kwargs)

        self.teacher_residuals = self.A @ self.teacher - self.b

        n = self.A.shape[0]
        self.noised_indices = noised_indices
        self.clean_indices = np.array([i for i in range(n) if i not in self.noised_indices])
        
        self.A_clean = self.A[self.clean_indices]
        self.min_sing_val_sq = np.linalg.svdvals(self.A_clean)[-1]**2

        self.clean_cardinality = self.A_clean.shape[0]
        self.tau = np.max(np.abs(self.teacher_residuals - self.true_residuals))
        self.max_residual = np.max(np.abs(self.true_residuals))

        self.alpha = self.min_sing_val_sq / self.clean_cardinality
        self.beta = 4 * self.tau
        self.gamma = 2 * self.max_residual * self.tau

        self.D = np.sqrt(self.beta**2 + 4 * self.alpha * self.gamma)
        self.r_plus = (self.beta + self.D) / (2 * self.alpha)
        self.r_minus = (self.beta - self.D) / (2 * self.alpha)

        self.get_bound_diff_eq_floor()
        self.get_bound_peter_paul_floor()
        
        self.fixed_shift = (self.beta + np.sqrt(self.beta**2 + 4 * self.alpha * self.gamma)) / (2 * self.alpha)
        self.get_bound_fixed_shift_floor()


        self._score()
        self.d0 = np.sqrt(self.parameter_errors[0]) - self.fixed_shift

        self.q = 1 - self.alpha + self.beta / (2 * self.fixed_shift)
        self.rho = max(np.sqrt(1 - self.alpha), self.q)

    def get_subset_indices(self, subset="all"):
        if subset == "all":
            return np.arange(self.n)
        elif subset == "clean":
            return self.clean_indices
        elif subset == "noised":
            return self.noised_indices
        else:
            raise NotImplementedError(f"Index subset {subset} not implemented for {self.name}")
    
    def select_index(self, indices):
        residuals = self.A[indices] @ self.x - self.b[indices]
        losses = (residuals - self.teacher_residuals[indices])**2 - self.teacher_residuals[indices]**2
        return indices[np.argmax(losses)]

    def get_bound(self):
        bound = (1 - self.min_sing_val_sq / self.A_clean.shape[0]) * np.sum((self.x - self.x_sol)**2)
        self.diagnostics["bounds"].append(bound)

        bound_fixed_shift = (self.fixed_shift + self.rho**self.k * self.d0)**2
        self.diagnostics["bound_fixed_shift"].append(bound_fixed_shift)

        test_val = self.rho**self.k * self.d0
        self.diagnostics["bound_fixed_shift_deviation"].append(test_val)

    def calc_diagnostics(self):
        self.get_bound()
        self.get_bound_diff_eq()

    def get_bound_diff_eq(self):
        error = self.parameter_errors[-1] # get last parameter error to bound min achievement time
        e0 = self.parameter_errors[0]
        def log_term(r):
            return r * np.log((np.sqrt(e0) - r) / (np.sqrt(error) - r))
        T = 2/self.D * (log_term(self.r_plus) - log_term(self.r_minus))
        self.diagnostics["bound_diff_eq"].append((T, error))

    def get_bound_diff_eq_floor(self):
        a = self.alpha
        b = self.beta
        g = self.gamma
        self.diagnostics["bound_diff_eq_floor"] = ((b + np.sqrt(b**2 + 4 * a * g)) / (2 * a))**2

    def get_bound_peter_paul_floor(self):
        a = self.alpha
        b = self.beta
        g = self.gamma
        self.diagnostics["bound_peter_paul_floor"] = 2 * (b**2 / a**2 + g / a)

    def get_bound_fixed_shift_floor(self):
        a = self.alpha
        b = self.beta
        g = self.gamma
        self.diagnostics["bound_fixed_shift_floor"] = self.fixed_shift**2

    def _score(self):
        residuals = (self.A @ self.x - self.b)**2
        
        self.scores.append(np.sum(residuals))
        self.parameter_errors.append(np.sum((self.x - self.x_sol)**2))
        self.diagnostics["scores_clean"].append(np.sum(residuals[self.clean_indices]))
        self.diagnostics["scores_noised"].append(np.sum(residuals[self.noised_indices]))


class Comparison:
    def __init__(self, solvers, name="comparison"):
        self.solvers = solvers
        self.plot_dir = "./plots"
        self.name = name

    def plot_parameter_errors(self, title=None, logx=False, logy=False):
        for sol in self.solvers:
            plt.plot(sol.parameter_errors, label=sol.name, lw=3)

        plt.xlabel("Iteration")
        plt.ylabel("Squared parameter error")
        plt.legend()
        plt.title("Parameter error" if title is None else title)

        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')

