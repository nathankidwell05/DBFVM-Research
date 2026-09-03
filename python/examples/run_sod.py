from pathlib import Path
import os

import jax.numpy as jnp
import numpy as np

from dbm2d.cases import make_grid, sod_initial_state
from dbm2d.config import SodConfig, SolverConfig
from dbm2d.equilibrium import feq_kt_d2v9
from dbm2d.exact_riemann import sod_exact
from dbm2d.moments import nonequilibrium_stress_heat_flux, primitive_from_f
from dbm2d.solver import run
from dbm2d.velocity_models import kt_d2v9


def main() -> None:
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))

    config = SolverConfig(nx=400, ny=8, lx=1.0, ly=0.02, tau=1.0e-4, t_end=0.15)
    sod = SodConfig(p_right=0.1)
    model = kt_d2v9()

    x, y, dx, dy = make_grid(config)
    _ = y
    f0, _ = sod_initial_state(config, sod, model)
    result = run(f0, config, model, dx, dy, jit=True)

    rho, ux, uy, temperature, pressure = primitive_from_f(result.f, model, config.gamma)
    feq = feq_kt_d2v9(rho, ux, uy, temperature, config.gamma, model)
    tne = nonequilibrium_stress_heat_flux(result.f, feq, model, ux, uy)

    mid = config.ny // 2
    x_np = np.asarray(x)
    exact = sod_exact(
        x_np,
        config.t_end,
        sod.diaphragm_x,
        (sod.rho_left, sod.u_left, sod.p_left),
        (sod.rho_right, sod.u_right, sod.p_right),
        config.gamma,
    )

    fields = {
        "x": x_np,
        "rho": np.asarray(rho[:, mid]),
        "ux": np.asarray(ux[:, mid]),
        "uy": np.asarray(uy[:, mid]),
        "pressure": np.asarray(pressure[:, mid]),
        "temperature": np.asarray(temperature[:, mid]),
        "rho_exact": exact["rho"],
        "ux_exact": exact["u"],
        "pressure_exact": exact["p"],
        "temperature_exact": exact["temperature"],
        "pi_xx": np.asarray(tne["pi_xx"][:, mid]),
        "pi_xy": np.asarray(tne["pi_xy"][:, mid]),
        "q_x": np.asarray(tne["q_x"][:, mid]),
    }
    np.savez(output_dir / "sod_centerline.npz", **fields)

    def l1(name: str, exact_name: str) -> float:
        return float(np.mean(np.abs(fields[name] - fields[exact_name])))

    y_variation = float(jnp.max(jnp.std(rho, axis=1)))
    summary = "\n".join(
        [
            "Finite-volume DBM Sod validation",
            f"model: {model.name}",
            f"grid: nx={config.nx}, ny={config.ny}",
            f"steps: {result.steps}",
            f"dt: {result.dt:.6e}",
            f"t_end: {result.time:.6e}",
            f"L1(rho): {l1('rho', 'rho_exact'):.6e}",
            f"L1(ux): {l1('ux', 'ux_exact'):.6e}",
            f"L1(pressure): {l1('pressure', 'pressure_exact'):.6e}",
            f"max y-row rho std: {y_variation:.6e}",
        ]
    )
    (output_dir / "sod_summary.txt").write_text(summary + "\n")
    print(summary)

    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
        plot_specs = [
            ("rho", "rho_exact", r"$\rho$"),
            ("ux", "ux_exact", r"$u_x$"),
            ("pressure", "pressure_exact", "p"),
            ("temperature", "temperature_exact", "T"),
        ]
        for ax, (name, exact_name, label) in zip(axes.ravel(), plot_specs):
            ax.plot(fields["x"], fields[name], label="DBM")
            ax.plot(fields["x"], fields[exact_name], "k--", label="exact")
            ax.set_xlabel("x")
            ax.set_ylabel(label)
            ax.grid(True, alpha=0.3)
            ax.legend()
        fig.savefig(output_dir / "sod_profiles.png", dpi=180)
    except Exception as exc:
        print(f"Skipping plot generation: {exc}")


if __name__ == "__main__":
    main()
