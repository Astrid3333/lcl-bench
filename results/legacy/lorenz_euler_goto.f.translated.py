import numpy as np

FLOAT_VARS = ['x', 'y', 'z']
STATE_NAMES = ['n', 'sigma', 'rho', 'beta', 'x', 'y', 'z', 'dt', 't']
LOOP_VAR = 'i'
MODE = 'goto_dispatch'

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
    _nb = 0
    _pc = 0
    while True:
        _nb += 1
        if _nb > 100000:
            raise RuntimeError("legacy step: goto budget exceeded (cyclic control flow)")
        if _pc == 0:
            dx = np.float64((sigma * (y - x)))
            dy = np.float64(((x * (rho - z)) - y))
            dz = np.float64(((x * y) - (beta * z)))
            x = np.float64((x + (dt * dx)))
            y = np.float64((y + (dt * dy)))
            z = np.float64((z + (dt * dz)))
            t = np.float64((t + dt))
            if ((abs(x) > np.float64(1e+30)) or (abs(y) > np.float64(1e+30))):
                _pc = 1
            else:
                _pc = 2
        elif _pc == 1:
            break
        elif _pc == 2:
            break
        else:
            break
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
