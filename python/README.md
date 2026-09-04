# Earlier JAX Prototype

This directory preserves the early JAX finite-volume D2V9 prototype created
before the project moved to the comment-heavy MATLAB workflow. It is retained
to document the project history and possible future vectorization ideas.

It is **not** the current reference solver, and its output is **not** used as
validation evidence in this repository. The validated work and reported
numerical results come from `matlab/d1v5/`.

To run the exploratory example:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python examples/run_sod.py
```

Generated files are written to `outputs/` and are not committed.
