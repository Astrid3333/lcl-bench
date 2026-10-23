"""Handles for named maps and ODEs (same math as octave/registry_map_step.m and registry_rhs.m)."""
import numpy as np


def map_step(name, p):
    if name == "logistic":
        r = p["r"]
        return lambda x: np.array([r * x[0] * (1.0 - x[0])])
    if name == "henon":
        a, b = p["a"], p["b"]
        return lambda x: np.array([1.0 - a * x[0] ** 2 + x[1], b * x[0]])
    if name == "clifford":
        a, b, c, d = p["a"], p["b"], p["c"], p["d"]
        return lambda x: np.array([
            np.sin(a * x[1]) + c * np.cos(a * x[0]),
            np.sin(b * x[0]) + d * np.cos(b * x[1]),
        ])
    if name == "tinkerbell":
        a, b, c, d = p["a"], p["b"], p["c"], p["d"]
        return lambda x: np.array([
            x[0] ** 2 - x[1] ** 2 + a * x[0] + b * x[1],
            2.0 * x[0] * x[1] + c * x[0] + d * x[1],
        ])
    raise KeyError(name)


def rhs(name, p):
    if name == "lorenz":
        s, r, b = p["sigma"], p["rho"], p["beta"]
        return lambda x: np.array([
            s * (x[1] - x[0]),
            x[0] * (r - x[2]) - x[1],
            x[0] * x[1] - b * x[2],
        ])
    if name == "rossler":
        a, b, c = p["a"], p["b"], p["c"]
        return lambda x: np.array([
            -x[1] - x[2],
            x[0] + a * x[1],
            b + x[2] * (x[0] - c),
        ])
    if name == "thomas":
        b = p["b"]
        return lambda x: np.array([
            np.sin(x[1]) - b * x[0],
            np.sin(x[2]) - b * x[1],
            np.sin(x[0]) - b * x[2],
        ])
    if name == "aizawa":
        a, b, c, d = p["a"], p["b"], p["c"], p["d"]
        e_, f_ = p["e"], p["f"]
        return lambda x: np.array([
            (x[2] - b) * x[0] - d * x[1],
            d * x[0] + (x[2] - b) * x[1],
            c + a * x[2] - x[2] ** 3 / 3.0
            - (x[0] ** 2 + x[1] ** 2) * (1.0 + e_ * x[2])
            + f_ * x[2] * x[0] ** 3,
        ])
    raise KeyError(name)
