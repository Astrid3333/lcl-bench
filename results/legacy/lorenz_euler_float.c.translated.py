import numpy as np

FLOAT_VARS = ['x', 'y', 'z']
STATE_NAMES = ['h', 'sigma', 'rho', 'beta', 'x', 'y', 'z']
LOOP_VAR = 'i'
MODE = 'structured'

def legacy_init():
    v = {}
    h = np.float32(np.float32(0.01))
    sigma = np.float32(np.float32(10.0))
    rho = np.float32(np.float32(28.0))
    beta = np.float32((np.float32(8.0) / np.float32(3.0)))
    x = np.float32(np.float32(1.0))
    y = np.float32(np.float32(1.0))
    z = np.float32(np.float32(1.0))
    v['h'] = h
    v['sigma'] = sigma
    v['rho'] = rho
    v['beta'] = beta
    v['x'] = x
    v['y'] = y
    v['z'] = z
    return v

def legacy_step(v, i=0):
    h = v['h']
    sigma = v['sigma']
    rho = v['rho']
    beta = v['beta']
    x = v['x']
    y = v['y']
    z = v['z']
    dx = np.float32((sigma * (y - x)))
    dy = np.float32(((x * (rho - z)) - y))
    dz = np.float32(((x * y) - (beta * z)))
    x = np.float32((x + (h * dx)))
    y = np.float32((y + (h * dy)))
    z = np.float32((z + (h * dz)))
    v['h'] = h
    v['sigma'] = sigma
    v['rho'] = rho
    v['beta'] = beta
    v['x'] = x
    v['y'] = y
    v['z'] = z
    return v

def legacy_steps(v=None):
    h = v['h']
    sigma = v['sigma']
    rho = v['rho']
    beta = v['beta']
    x = v['x']
    y = v['y']
    z = v['z']
    _s = (0 + 1)
    _e = 50000
    return max(0, int(_e) - int(_s) + 1)
