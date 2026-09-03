# Earlier JAX prototype

This is the original 2D finite-volume DBM prototype. It uses a KT-style D2V9
model on a planar Sod problem and includes an exact Euler Riemann comparison.
It is useful for architecture and regression work, but the MATLAB D1V5 solver
is the current reference implementation.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
python examples/run_sod.py
pytest
```

Generated files are written to `outputs/` and are not committed.
