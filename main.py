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

        self.rng = np.random.default_rng(seed=seed)

        self.scores = []
        self._score()
        self.bounds = [float('nan')]

        self.diagnostics = defaultdict(list)

        self.name = type(self).__name__

    def reset(self):
        self.x = self.x0.copy()
        self.k = 0
        self.scores = []
        self._score()

        self.diagnostics = defaultdict(list)

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

    def calc_diagnostics(self):
        return None
    
    def solve(self):
        while not self.converged():
            print(f'Iterating (k = {self.k})', end='\r')
            self.calc_diagnostics()
            self._step()
            self._score()
        
        return self.x

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


class MillerOpt(KaczmarzSolver):
    def __init__(self, *args, teacher, noised_indices=[], bound_type=0, **kwargs):
        self.teacher = teacher

        super().__init__(*args, **kwargs)

        self.teacher_residuals = self.A @ self.teacher - self.b
        self.true_residuals = self.A @ self.x_sol - self.b

        n = self.A.shape[0]
        self.A_clean = self.A[[i for i in range(n) if i not in noised_indices]]
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

        self.get_time_bound_floor()
    
    def select_index(self, indices):
        residuals = self.A[indices] @ self.x - self.b[indices]
        losses = (residuals - self.teacher_residuals[indices])**2 - self.teacher_residuals[indices]**2
        return indices[np.argmax(losses)]

    def get_bound(self):
        bound = (1 - self.min_sing_val_sq / self.A_clean.shape[0]) * np.sum((self.x - self.x_sol)**2)
        self.diagnostics["bounds"].append(bound)

    def calc_diagnostics(self):
        self.get_bound()
        self.get_time_bound()

    def get_time_bound(self):
        error = self.scores[-1] # get last score to bound min achievement time
        e0 = self.scores[0]
        def log_term(r):
            return r * np.log((np.sqrt(e0) - r) / (np.sqrt(error) - r))
        T = 2/self.D * (log_term(self.r_plus) - log_term(self.r_minus))
        self.diagnostics["time bounds"].append((T, error))

    def get_time_bound_floor(self):
        a = self.alpha
        b = self.beta
        g = self.gamma
        self.tb_floor = ((b + np.sqrt(b**2 + 4 * a * g)) / (2 * a))**2
        return self.tb_floor


class Comparison:
    def __init__(self, solvers, name="comparison"):
        self.solvers = solvers
        self.plot_dir = "./plots"
        self.name = name

    def plot_residuals(self, title=None, name_suffix='', logx=False, logy=False):
        for sol in self.solvers:
            plt.plot(sol.scores, label=sol.name)
            if isinstance(sol, MillerOpt):
                plt.axhline(np.sum(sol.teacher_residuals**2), label=f"{sol.name} teacher", ls="--")
            

        plt.xlabel("Iteration")
        plt.ylabel("Sum of Squared Errors (SSE)")
        plt.legend()
        plt.title("Comparing convergence of Kaczmarz methods" if title is None else title)

        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')

        filename = self.name + name_suffix + ".png"
        plt.savefig(os.path.join(self.plot_dir, filename))

    # def plot_time_bounds(self, title=None, name_suffix='', logx=False, logy=False):
    #     for sol in self.solvers:
    #         plt.plot(sol.time_bounds, label=sol.name)
    #         if isinstance(sol, MillerOpt):
    #             plt.axhline(np.sum(sol.teacher_residuals**2), label=f"{sol.name} teacher", ls="--")
        