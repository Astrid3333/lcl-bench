import numpy as np

FLOAT_VARS = ['x', 'y', 'z']
STATE_NAMES = ['n', 'sigma', 'rho', 'beta', 'x', 'y', 'z', 'dt', 't']
LOOP_VAR = 'i'

def legacy_init():
    v = {}
    n = int(50000)
    sigma = np.float32(np.float32(10.0))
    rho = np.float32(np.float32(28.0))
    beta = np.float32(np.float32(2.6666667))
    x = np.float32(np.float32(1.0))
    y = np.float32(np.float32(1.0))
    z = np.float32(np.float32(1.0))
    dt = np.float32(np.float32(0.01))
    t = np.float32(np.float32(0.0))
    v['n'] = n
    v['sigma'] = sigma
    v['rho'] = rho
    v['beta'] = beta
    v['x'] = x
    v['y'] = y
    v['z'] = z
    v['dt'] = dt
    v['t'] = t
    return v

def legacy_step(v, i=0):
    n = v['n']
    sigma = v['sigma']
    rho = v['rho']
    beta = v['beta']
    x = v['x']
    y = v['y']
    z = v['z']
    dt = v['dt']
    t = v['t']
    dx = np.float32((sigma * (y - x)))
    dy = np.float32(((x * (rho - z)) - y))
    dz = np.float32(((x * y) - (beta * z)))
    x = np.float32((x + (dt * dx)))
    y = np.float32((y + (dt * dy)))
    z = np.float32((z + (dt * dz)))
    t = np.float32((t + dt))
    v['n'] = n
    v['sigma'] = sigma
    v['rho'] = rho
    v['beta'] = beta
    v['x'] = x
    v['y'] = y
    v['z'] = z
    v['dt'] = dt
    v['t'] = t
    return v

def legacy_steps(v=None):
    n = v['n']
    sigma = v['sigma']
    rho = v['rho']
    beta = v['beta']
    x = v['x']
    y = v['y']
    z = v['z']
    dt = v['dt']
    t = v['t']
    _s = 1
    _e = n
    return max(0, int(_e) - int(_s) + 1)

def legacy_final(v):
    x = v['x']
    y = v['y']
    z = v['z']
    return [x, y, z]
