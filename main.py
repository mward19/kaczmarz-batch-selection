from typing import List
import os

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

        self.name = type(self).__name__

    def reset(self):
        self.x = self.x0.copy()
        self.k = 0
        self.scores = []
        self._score()
        self.bounds = [float('nan')]

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
        score = np.mean((self.x - self.x_sol)**2)
        self.scores.append(score)

    def get_bound(self):
        raise NotImplementedError()
    
    def _calc_bound(self):
        self.bounds.append(self.get_bound())
    
    def solve(self):
        while not self.converged():
            print(f'Iterating (k = {self.k})', end='\r')
            self._calc_bound()
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

    def get_bound(self):
        # Expected convergence rate for randomized Kaczmarz (rows are normalized,
        # so the squared Frobenius norm is just the row count)
        return (1 - self.min_sing_val_sq / self.A.shape[0]) * np.sum((self.x - self.x_sol)**2)


class Motzkins(KaczmarzSolver):
    def select_index(self, indices):
        residuals = self.A[indices] @ self.x - self.b[indices]
        return indices[np.argmax(residuals**2)]
    
    def get_bound(self):
        # For the upcoming iteration. Lemma 1
        max_residual = np.max(np.abs(self.A @ self.x - self.b))
        if max_residual > 4 * self.max_e: # not close to sol'n
            return np.sum((self.x - self.x_sol)**2) - 1/2 * max_residual**2
        else:
            return float('nan') # temporary placeholder


class MillerOpt(KaczmarzSolver):
    def __init__(self, *args, teacher, noised_indices=[], **kwargs):
        self.teacher = teacher

        super().__init__(*args, **kwargs)

        self.teacher_residuals = (self.A @ self.teacher - self.b)

        n = self.A.shape[0]
        self.A_clean = self.A[[i for i in range(n) if i not in noised_indices]]
        self.min_sing_val_sq = np.linalg.svdvals(self.A_clean)[-1]**2
    
    def select_index(self, indices):
        residuals = self.A[indices] @ self.x - self.b[indices]
        losses = (residuals - self.teacher_residuals[indices])**2 - self.teacher_residuals[indices]**2
        return indices[np.argmax(losses)]

    def get_bound(self):
        bound = (1 - self.min_sing_val_sq / self.A_clean.shape[0]) * np.sum((self.x - self.x_sol)**2)
        return bound


class Comparison:
    def __init__(self, solvers, name="comparison"):
        self.solvers = solvers
        self.plot_dir = "./plots"
        self.name = name

    def plot_residuals(self, title=None, name_suffix='', logx=False, logy=False):
        for sol in self.solvers:
            plt.plot(sol.scores, label=sol.name)
            if isinstance(sol, MillerOpt):
                plt.axhline(np.mean(sol.teacher_residuals**2), label=f"{sol.name} teacher", ls="--")
            

        plt.xlabel("Iteration")
        plt.ylabel("Mean Squared Error (MSE)")
        plt.legend()
        plt.title("Comparing convergence of Kaczmarz methods" if title is None else title)

        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')

        filename = self.name + name_suffix + ".png"
        plt.savefig(os.path.join(self.plot_dir, filename))


if __name__ == "__main__":
    # Generate system
    rng = np.random.default_rng(seed=1)
    n, d = (100, 50)
    A, x, b = get_system(rng, n, d)

    x_teacher = x + 0.5 * rng.normal(size=d)

    x0 = rng.normal(size=d)

    common_params = dict(
        x0=x0,
        A=A,
        b=b,
        x_sol=x,
        batch_size=10,
        maxiter=5000
    )
    
    sols = [
        MillerOpt(**common_params, teacher=x_teacher),
        RandomizedKaczmarz(**common_params),
    ]

    for sol in sols:
        sol.solve()

    comparison = Comparison(sols)
    comparison.plot_residuals(
        title=f"Comparing Kaczmarz convergence (batch size {common_params["batch_size"]})",
        logx=True,
        logy=False
    )
