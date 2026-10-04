“”“Helper functions for the yam peel adsorption analysis.

All models are defined in their non-linear (untransformed) form. Linearised
forms appear only where a notebook explicitly compares them with the
non-linear fit.
“””

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

R_GAS = 8.314  # J mol-1 K-1

# –––––––––––––––––––––––––––––––– data

def load_readings(path=“data/aas_readings.csv”):
“”“Load the raw AAS readings exactly as reported by the laboratory.”””
return pd.read_csv(path)

def concentration(reading, dilution_factor):
“”“Filtrate concentration (mg/L) on the reading x DF convention.

```
The laboratory's working notes were not available, so this convention is
inferred from the report layout. It affects absolute concentrations and
uptake values, but not removal percentages calculated against a control
analysed at the same dilution.
"""
return np.asarray(reading, dtype=float) * np.asarray(dilution_factor, dtype=float)
```

def removal_pct(c0, ce):
“”“Percentage removal from initial and residual concentrations.”””
c0 = np.asarray(c0, dtype=float)
ce = np.asarray(ce, dtype=float)
return (c0 - ce) / c0 * 100.0

def uptake(c0, ce, volume_L, mass_g):
“”“Amount adsorbed per gram, qe = (C0 - Ce) * V / m, in mg/g.”””
c0 = np.asarray(c0, dtype=float)
ce = np.asarray(ce, dtype=float)
volume_L = np.asarray(volume_L, dtype=float)
return (c0 - ce) * volume_L / np.asarray(mass_g, dtype=float)

# ———————————————————— isotherms

def langmuir(ce, qmax, kl):
return qmax * kl * ce / (1.0 + kl * ce)

def freundlich(ce, kf, inv_n):
return kf * np.power(ce, inv_n)

def temkin(ce, b, at):
return b * np.log(at * ce)

# ———————————————————–– kinetics

def plateau(t, qe):
“”“Constant uptake: the system has already equilibrated before the first
measured time. One parameter, no time dependence.”””
return np.full_like(np.asarray(t, dtype=float), qe)

def pfo(t, qe, k1):
return qe * (1.0 - np.exp(-k1 * np.asarray(t, dtype=float)))

def pso(t, qe, k2):
t = np.asarray(t, dtype=float)
return k2 * qe ** 2 * t / (1.0 + k2 * qe * t)

def elovich(t, alpha, beta):
t = np.asarray(t, dtype=float)
return (1.0 / beta) * np.log(1.0 + alpha * beta * t)

def intraparticle(t, kid, c):
return kid * np.sqrt(np.asarray(t, dtype=float)) + c

# —————————————————— fit diagnostics

def aicc(sse, n, k):
“”“Small-sample Akaike information criterion.

```
Lower is better. A difference of more than about 2 units is usually
treated as meaningful support for the lower model.
"""
if n - k - 1 <= 0:
    return np.nan
return n * np.log(sse / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
```

def fit_model(func, x, y, p0, bounds=(-np.inf, np.inf), maxfev=200000):
“”“Least-squares fit returning parameters and fit diagnostics.

```
Returns a dict. If the optimiser fails to converge, 'converged' is False
and the parameters are None. Non-convergence is reported, not hidden.
"""
x = np.asarray(x, dtype=float)
y = np.asarray(y, dtype=float)
n = len(y)
k = len(p0)
try:
    popt, _ = curve_fit(func, x, y, p0=p0, bounds=bounds, maxfev=maxfev)
except Exception as exc:
    return {"converged": False, "message": str(exc), "params": None,
            "r2": np.nan, "rmse": np.nan, "aicc": np.nan, "sse": np.nan}
pred = func(x, *popt)
resid = y - pred
sse = float(np.sum(resid ** 2))
sst = float(np.sum((y - y.mean()) ** 2))
r2 = 1.0 - sse / sst if sst > 0 else np.nan
return {"converged": True, "message": "", "params": popt,
        "r2": r2, "rmse": float(np.sqrt(sse / n)),
        "aicc": aicc(sse, n, k), "sse": sse}
```

def fit_plateau(y):
“”“The plateau fit is the mean of the observations (one parameter).”””
y = np.asarray(y, dtype=float)
n = len(y)
qe = float(y.mean())
sse = float(np.sum((y - qe) ** 2))
sst = float(np.sum((y - y.mean()) ** 2))
return {“converged”: True, “message”: “”, “params”: np.array([qe]),
“r2”: 0.0 if sst == 0 else 1.0 - sse / sst,
“rmse”: float(np.sqrt(sse / n)),
“aicc”: aicc(sse, n, 1), “sse”: sse}
