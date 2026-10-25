import numpy as np

FLOAT_VARS = ['x', 'y', 'z']
STATE_NAMES = ['n', 'sigma', 'rho', 'beta', 'x', 'y', 'z', 'dt', 't']
LOOP_VAR = 'i'

def legacy_init():
    v = {}
    n = int(100000)
    sigma = np.float64(np.float64(10.0))
    rho = np.float64(np.float64(28.0))
    beta = np.float64(np.float64(2.666666666666667))
    x = np.float64(np.float64(1.0))
    y = np.float64(np.float64(1.0))
    z = np.float64(np.float64(1.0))
    dt = np.float64(np.float64(0.01))
    t = np.float64(np.float64(0.0))
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
    k1x = np.float64(((dt * sigma) * (y - x)))
    k1y = np.float64((dt * ((x * (rho - z)) - y)))
    k1z = np.float64((dt * ((x * y) - (beta * z))))
    k2x = np.float64(((dt * sigma) * ((y + (np.float64(0.5) * k1y)) - (x + (np.float64(0.5) * k1x)))))
    k2y = np.float64((dt * (((x + (np.float64(0.5) * k1x)) * (rho - (z + (np.float64(0.5) * k1z)))) - (y + (np.float64(0.5) * k1y)))))
    k2z = np.float64((dt * (((x + (np.float64(0.5) * k1x)) * (y + (np.float64(0.5) * k1y))) - (beta * (z + (np.float64(0.5) * k1z))))))
    k3x = np.float64(((dt * sigma) * ((y + (np.float64(0.5) * k2y)) - (x + (np.float64(0.5) * k2x)))))
    k3y = np.float64((dt * (((x + (np.float64(0.5) * k2x)) * (rho - (z + (np.float64(0.5) * k2z)))) - (y + (np.float64(0.5) * k2y)))))
    k3z = np.float64((dt * (((x + (np.float64(0.5) * k2x)) * (y + (np.float64(0.5) * k2y))) - (beta * (z + (np.float64(0.5) * k2z))))))
    k4x = np.float64(((dt * sigma) * ((y + k3y) - (x + k3x))))
    k4y = np.float64((dt * (((x + k3x) * (rho - (z + k3z))) - (y + k3y))))
    k4z = np.float64((dt * (((x + k3x) * (y + k3y)) - (beta * (z + k3z)))))
    x = np.float64((x + ((((k1x + (np.float64(2.0) * k2x)) + (np.float64(2.0) * k3x)) + k4x) / np.float64(6.0))))
    y = np.float64((y + ((((k1y + (np.float64(2.0) * k2y)) + (np.float64(2.0) * k3y)) + k4y) / np.float64(6.0))))
    z = np.float64((z + ((((k1z + (np.float64(2.0) * k2z)) + (np.float64(2.0) * k3z)) + k4z) / np.float64(6.0))))
    t = np.float64((t + dt))
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
