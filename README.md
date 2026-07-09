# ai-rd-assignment
Recovering hidden parameters (θ, M, X) of a parametric curve from unordered sample points using iterative orthogonal distance regression — R&amp;D/AI assignment
# AI R&D Assignment — Parametric Curve Fit

## Problem
Recover `theta`, `M`, `X` in:

```
x(t) = t*cos(theta) - e^(M|t|)*sin(0.3t)*sin(theta) + X
y(t) = 42 + t*sin(theta) + e^(M|t|)*sin(0.3t)*cos(theta)
```

for `6 < t < 60`, `0° < theta < 50°`, `-0.05 < M < 0.05`, `0 < X < 100`, given
1500 `(x, y)` points sampled from the curve (`xy_data.csv`).

## Approach

**1. Inspect the data.**
A scatter plot of the raw `(x, y)` pairs traces a single clean, smooth curve
with no scatter/noise — so the points are exact (or near-exact) samples of
the true curve, just in shuffled row order (the CSV order does **not**
correspond to increasing `t`). That ruled out simply assuming row `i` maps to
`t_i = linspace(6, 60, 1500)[i]`.

**2. Treat it as curve fitting to an unordered point cloud.**
Since each point's own `t` is unknown, this is an orthogonal-distance
regression problem: find `(theta, M, X)` that minimizes the total distance
from each data point to *its closest point* on the model curve, over all
valid `t` in `[6, 60]`.

**3. Alternating optimization (ICP-style), iterated to convergence:**
- **E-step:** for the current `(theta, M, X)`, find the best-matching `t`
  for every data point (coarse search over a dense `t` grid, then refined
  with a bounded 1-D local minimizer).
- **M-step:** with those per-point `t` values fixed, refit `(theta, M, X)`
  with nonlinear least squares (`scipy.optimize.least_squares`,
  Levenberg–Marquardt/trust-region, respecting the given parameter bounds).
- Repeat. Each round only tightens the fit, and the cost (sum of squared
  residuals) drops monotonically from ~7300 to effectively 0.

**4. Convergence.**
After ~35–40 iterations the parameters lock onto clean, "round" values and
the max per-point residual falls to numerical noise (~1e-5), meaning the
model curve passes through every data point almost exactly:

| Parameter | Recovered value |
|---|---|
| theta | **30°** (0.5235987756 rad) |
| M | **0.03** |
| X | **55** |

**5. Verification.**
Plotting the fitted curve `(x(t), y(t))` for `t ∈ [6, 60]` on top of the
1500 data points shows a perfect overlay (`verify_fit.png`).

## Files
- `fit.py` — full fitting pipeline (run with `xy_data.csv` in the same folder)
- `xy_data.csv` — provided dataset
- `verify_fit.png` — fitted curve overlaid on the given data points

## Result — required LaTeX / Desmos submission

```
\left(t*\cos(0.5235987756)-e^{0.03\left|t\right|}\cdot\sin(0.3t)\sin(0.5235987756)+55,42+t*\sin(0.5235987756)+e^{0.03\left|t\right|}\cdot\sin(0.3t)\cos(0.5235987756)\right)
```
with domain `6 ≤ t ≤ 60`.

Equivalently:
- **theta = 30°**
- **M = 0.03**
- **X = 55**
