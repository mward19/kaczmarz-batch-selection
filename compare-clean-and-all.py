import marimo

__generated_with = "0.25.2"
app = marimo.App()


@app.cell
def _():
    from main import TeacherOpt, Comparison, RandomizedKaczmarz, get_system
    from matplotlib import pyplot as plt
    import marimo as mo

    return Comparison, TeacherOpt, get_system, mo, plt


@app.cell
def _(mo):
    import numpy as np

    n = 1000
    d = 50
    prop_noised = 0.1
    b_noise_scale = 0.2
    teacher_noise_scale = 0.1

    iters = 5000
    seed = 1

    _steps = np.round(np.geomspace(0.01, 300, 1000), decimals=2)
    slider = mo.ui.slider(
        steps=_steps, 
        label="Variance of 0-mean initial guess", 
        value=_steps[-1],
        full_width=True
    )
    return (
        b_noise_scale,
        d,
        iters,
        n,
        np,
        prop_noised,
        seed,
        slider,
        teacher_noise_scale,
    )


@app.cell
def _(mo, slider):
    mo.hstack([slider, slider.value])
    return


@app.cell
def _(
    TeacherOpt,
    b_noise_scale,
    d,
    get_system,
    iters,
    n,
    np,
    prop_noised,
    seed,
    slider,
    teacher_noise_scale,
):
    x0_scale = slider.value

    # Generate system
    rng = np.random.default_rng(seed=seed)
    n_noised = int(n * prop_noised)
    A, x, b = get_system(rng, n, d)

    noised_indices = rng.choice(n, size=n_noised, replace=False)
    b_noised = b
    b_noised[noised_indices] += rng.normal(scale=b_noise_scale, size=n_noised)

    x_teacher = x + rng.normal(size=d, scale=teacher_noise_scale)

    x0 = rng.normal(scale=x0_scale**0.5, size=d)

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
    return common_params, sols


@app.cell
def _(Comparison, b_noise_scale, common_params, iters, np, plt, sols):
    _comparison = Comparison(sols)
    _comparison.plot_parameter_errors(
        title=f"Comparing Kaczmarz convergence\n(batch size {common_params['batch_size']}, noise std {b_noise_scale})",
        logx=True,
        logy=True
    )

    _bound_diff_eq = np.array(sols[0].diagnostics["bound_diff_eq"])

    plt.plot(_bound_diff_eq[:, 0], _bound_diff_eq[:, 1], color="tab:green", label="Bound (DE)", lw=1)
    plt.axhline(sols[0].diagnostics["bound_diff_eq_floor"], label="Floor (DE)", ls="--", color="tab:green", lw=1)

    # plt.plot(sols[0].diagnostics["bounds"], color="tab:red", label="Bound (PP)")
    _iter_range = np.arange(iters)
    _e0 = sols[0].parameter_errors[0]
    _bound_peter_paul_floor = sols[0].diagnostics["bound_peter_paul_floor"]
    plt.plot(_iter_range + 1, (1 - sols[0].alpha / 2)**_iter_range * (_e0 - _bound_peter_paul_floor) + _bound_peter_paul_floor, color="tab:red", label="Bound (PP)", lw=1)
    # plt.axhline(_bound_peter_paul_floor, color="tab:red", ls="--", label="Floor (PP)", lw=1)

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


@app.cell
def _(Comparison, b_noise_scale, common_params, iters, np, plt, sols):
    _comparison = Comparison(sols)
    _comparison.plot_parameter_errors(
        title=f"Parameter error\n(batch size {common_params['batch_size']}, noise std {b_noise_scale})",
        logx=True,
        logy=True
    )

    _bound_diff_eq = np.array(sols[0].diagnostics["bound_diff_eq"])

    plt.plot(_bound_diff_eq[:, 0], _bound_diff_eq[:, 1], color="tab:green", label="Bound (DE)", lw=1)
    plt.axhline(sols[0].diagnostics["bound_diff_eq_floor"], label="Floor (DE)", ls="--", color="tab:green", lw=1)

    # plt.plot(sols[0].diagnostics["bounds"], color="tab:red", label="Bound (PP)")
    iter_range = np.arange(iters)
    _e0 = sols[0].parameter_errors[0]
    _bound_peter_paul_floor = sols[0].diagnostics["bound_peter_paul_floor"]
    plt.plot(iter_range + 1, (1 - sols[0].alpha / 2)**iter_range * (_e0 - _bound_peter_paul_floor) + _bound_peter_paul_floor, color="tab:red", label="Bound (PP)", lw=1)
    plt.axhline(_bound_peter_paul_floor, color="tab:red", ls="--", label="Floor (PP)", lw=1)

    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.show()

    # print(f"alpha: {sols[0].alpha}")
    # print(f"beta:  {sols[0].beta}")
    # print(f"gamma: {sols[0].gamma}")
    return


@app.cell
def _(sols):
    print(f"final res all:   {sols[0].scores[-1]}")
    print(f"final res clean: {sols[0].diagnostics["scores_clean"][-1]}")
    return


@app.cell
def _(iters, np, plt, sols):
    _iter_range = np.arange(iters)
    plt.semilogy(_iter_range, sols[0].parameter_errors[:iters], label="Squared parameter error", ls="-", color="tab:purple", lw=2)
    return


@app.cell
def _(Comparison, b_noise_scale, common_params, iters, np, plt, sols):
    _comparison = Comparison(sols)
    _comparison.plot_parameter_errors(
        title=f"Comparing Kaczmarz convergence\n(batch size {common_params['batch_size']}, noise std {b_noise_scale})",
        logx=True,
        logy=True
    )

    _bound_diff_eq = np.array(sols[0].diagnostics["bound_diff_eq"])

    plt.plot(_bound_diff_eq[:, 0], (_bound_diff_eq[:, 1])**0.5, color="tab:green", label="Bound (DE)", lw=1)
    plt.axhline((sols[0].diagnostics["bound_diff_eq_floor"])**0.5, label="Floor (DE)", ls="--", color="tab:green", lw=1)

    # plt.plot(sols[0].diagnostics["bounds"], color="tab:red", label="Bound (PP)")
    _iter_range = np.arange(iters)
    _e0 = sols[0].parameter_errors[0]
    _bound_peter_paul_floor = sols[0].diagnostics["bound_peter_paul_floor"]
    plt.plot(_iter_range + 1, ((1 - sols[0].alpha / 2)**_iter_range * (_e0 - _bound_peter_paul_floor) + _bound_peter_paul_floor)**0.5, color="tab:red", label="Bound (PP)", lw=1)
    # plt.axhline(_bound_peter_paul_floor, color="tab:red", ls="--", label="Floor (PP)", lw=1)

    _bound_fixed_shift_floor = sols[0].diagnostics["bound_fixed_shift_floor"]
    plt.axhline(_bound_fixed_shift_floor**0.5, color="tab:purple", ls="--", label="Floor (fixed point shift)", lw=1)

    plt.plot(_iter_range, np.array(sols[0].diagnostics["bound_fixed_shift"])**0.5, label="Bound (fixed point shift)", ls="-", color="tab:purple", lw=2)

    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))

    plt.show()

    # print(f"alpha: {sols[0].alpha}")
    # print(f"beta:  {sols[0].beta}")
    # print(f"gamma: {sols[0].gamma}")
    return


if __name__ == "__main__":
    app.run()
