import marimo

__generated_with = "0.25.2"
app = marimo.App()


@app.cell
def _():
    import numpy as np
    from main import TeacherOpt, Comparison, RandomizedKaczmarz, get_system, get_worst_system
    from matplotlib import pyplot as plt
    import marimo as mo

    return Comparison, TeacherOpt, get_worst_system, np, plt


@app.cell
def _(TeacherOpt, get_worst_system, np):
    # Generate system
    rng = np.random.default_rng(seed=0)
    iters = 5
    n_noised = 1
    A, x, b, x_teacher = get_worst_system()

    print(x)
    n, d = A.shape

    noised_indices = np.array([4])
    b_noised = b
    # b_noised[noised_indices] += rng.normal(scale=b_noise_scale, size=n_noised)

    R = 1000
    x0 = np.full(2, R).astype(float)

    common_params = dict(
        x0=x0,
        A=A,
        b=b_noised,
        x_sol=x,
        batch_size=n,
        maxiter=iters
    )

    sols = [
        TeacherOpt(**common_params, teacher=x_teacher, noised_indices=noised_indices),
        # RandomizedKaczmarz(**common_params),
    ]

    for sol in sols:
        sol.solve()
    return iters, sols


@app.cell
def _(Comparison, iters, np, plt, sols):
    _comparison = Comparison(sols)
    _comparison.plot_parameter_errors(
        title=f"Comparing Kaczmarz convergence",
        logx=False,
        logy=True
    )

    # _bound_diff_eq = np.array(sols[0].diagnostics["bound_diff_eq"])

    # plt.plot(_bound_diff_eq[:, 0], _bound_diff_eq[:, 1], color="tab:green", label="Bound (DE)", lw=1)
    # plt.axhline(sols[0].diagnostics["bound_diff_eq_floor"], label="Floor (DE)", ls="--", color="tab:green", lw=1)

    # # plt.plot(sols[0].diagnostics["bounds"], color="tab:red", label="Bound (PP)")
    _iter_range = np.arange(iters)
    _e0 = sols[0].parameter_errors[0]
    # _bound_peter_paul_floor = sols[0].diagnostics["bound_peter_paul_floor"]
    # plt.plot(_iter_range + 1, (1 - sols[0].alpha / 2)**_iter_range * (_e0 - _bound_peter_paul_floor) + _bound_peter_paul_floor, color="tab:red", label="Bound (PP)", lw=1)
    # # plt.axhline(_bound_peter_paul_floor, color="tab:red", ls="--", label="Floor (PP)", lw=1)

    _bound_fixed_shift_floor = sols[0].diagnostics["bound_fixed_shift_floor"]
    plt.axhline(_bound_fixed_shift_floor, color="tab:purple", ls="--", label="Floor (fixed point shift)", lw=1)

    plt.plot(_iter_range, sols[0].diagnostics["bound_fixed_shift"], label="Bound (fixed point shift)", ls="-", color="tab:purple", lw=2)

    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))

    plt.show()

    # print(f"alpha: {sols[0].alpha}")
    # print(f"beta:  {sols[0].beta}")
    # print(f"gamma: {sols[0].gamma}")
    return


@app.cell
def _(np, sols):
    print(f"Rho: {sols[0].rho}")
    print(f"q: {sols[0].q}")
    print(f"sqrt(1 - alpha): {np.sqrt(1 - sols[0].alpha)}")
    return


@app.cell
def _(sols):
    sols[0].rho
    return


if __name__ == "__main__":
    app.run()
