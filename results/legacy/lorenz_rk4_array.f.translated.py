import numpy as np

FLOAT_VARS = ['y[0]', 'y[1]', 'y[2]']
STATE_NAMES = ['n', 'sigma', 'rho', 'beta', 'dt']
LOOP_VAR = 'i'
MODE = 'structured'

def legacy_init():
    v = {}
    y = [np.float64(0.0)] * 3
    k1 = [np.float64(0.0)] * 3
    k2 = [np.float64(0.0)] * 3
    k3 = [np.float64(0.0)] * 3
    k4 = [np.float64(0.0)] * 3
    n = int(100000)
    sigma = np.float64(np.float64(10.0))
    rho = np.float64(np.float64(28.0))
    beta = np.float64(np.float64(2.666666666666667))
    y[0] = np.float64(np.float64(1.0))
    y[1] = np.float64(np.float64(1.0))
    y[2] = np.float64(np.float64(1.0))
    dt = np.float64(np.float64(0.01))
    v['n'] = n
    v['sigma'] = sigma
    v['rho'] = rho
    v['beta'] = beta
    v['dt'] = dt
    v['y[0]'] = y[0]
    v['y[1]'] = y[1]
    v['y[2]'] = y[2]
    v['k1[0]'] = k1[0]
    v['k1[1]'] = k1[1]
    v['k1[2]'] = k1[2]
    v['k2[0]'] = k2[0]
    v['k2[1]'] = k2[1]
    v['k2[2]'] = k2[2]
    v['k3[0]'] = k3[0]
    v['k3[1]'] = k3[1]
    v['k3[2]'] = k3[2]
    v['k4[0]'] = k4[0]
    v['k4[1]'] = k4[1]
    v['k4[2]'] = k4[2]
    return v

def legacy_step(v, i=0):
    n = v['n']
    sigma = v['sigma']
    rho = v['rho']
    beta = v['beta']
    dt = v['dt']
    y = [v['y[0]'], v['y[1]'], v['y[2]']]
    k1 = [v['k1[0]'], v['k1[1]'], v['k1[2]']]
    k2 = [v['k2[0]'], v['k2[1]'], v['k2[2]']]
    k3 = [v['k3[0]'], v['k3[1]'], v['k3[2]']]
    k4 = [v['k4[0]'], v['k4[1]'], v['k4[2]']]
    k1[0] = np.float64(((dt * sigma) * (y[1] - y[0])))
    k1[1] = np.float64((dt * ((y[0] * (rho - y[2])) - y[1])))
    k1[2] = np.float64((dt * ((y[0] * y[1]) - (beta * y[2]))))
    k2[0] = np.float64(((dt * sigma) * ((y[1] + (np.float64(0.5) * k1[1])) - (y[0] + (np.float64(0.5) * k1[0])))))
    k2[1] = np.float64((dt * (((y[0] + (np.float64(0.5) * k1[0])) * (rho - (y[2] + (np.float64(0.5) * k1[2])))) - (y[1] + (np.float64(0.5) * k1[1])))))
    k2[2] = np.float64((dt * (((y[0] + (np.float64(0.5) * k1[0])) * (y[1] + (np.float64(0.5) * k1[1]))) - (beta * (y[2] + (np.float64(0.5) * k1[2]))))))
    k3[0] = np.float64(((dt * sigma) * ((y[1] + (np.float64(0.5) * k2[1])) - (y[0] + (np.float64(0.5) * k2[0])))))
    k3[1] = np.float64((dt * (((y[0] + (np.float64(0.5) * k2[0])) * (rho - (y[2] + (np.float64(0.5) * k2[2])))) - (y[1] + (np.float64(0.5) * k2[1])))))
    k3[2] = np.float64((dt * (((y[0] + (np.float64(0.5) * k2[0])) * (y[1] + (np.float64(0.5) * k2[1]))) - (beta * (y[2] + (np.float64(0.5) * k2[2]))))))
    k4[0] = np.float64(((dt * sigma) * ((y[1] + k3[1]) - (y[0] + k3[0]))))
    k4[1] = np.float64((dt * (((y[0] + k3[0]) * (rho - (y[2] + k3[2]))) - (y[1] + k3[1]))))
    k4[2] = np.float64((dt * (((y[0] + k3[0]) * (y[1] + k3[1])) - (beta * (y[2] + k3[2])))))
    y[0] = np.float64((y[0] + ((((k1[0] + (np.float64(2.0) * k2[0])) + (np.float64(2.0) * k3[0])) + k4[0]) / np.float64(6.0))))
    y[1] = np.float64((y[1] + ((((k1[1] + (np.float64(2.0) * k2[1])) + (np.float64(2.0) * k3[1])) + k4[1]) / np.float64(6.0))))
    y[2] = np.float64((y[2] + ((((k1[2] + (np.float64(2.0) * k2[2])) + (np.float64(2.0) * k3[2])) + k4[2]) / np.float64(6.0))))
    v['n'] = n
    v['sigma'] = sigma
    v['rho'] = rho
    v['beta'] = beta
    v['dt'] = dt
    v['y[0]'] = y[0]
    v['y[1]'] = y[1]
    v['y[2]'] = y[2]
    v['k1[0]'] = k1[0]
    v['k1[1]'] = k1[1]
    v['k1[2]'] = k1[2]
    v['k2[0]'] = k2[0]
    v['k2[1]'] = k2[1]
    v['k2[2]'] = k2[2]
    v['k3[0]'] = k3[0]
    v['k3[1]'] = k3[1]
    v['k3[2]'] = k3[2]
    v['k4[0]'] = k4[0]
    v['k4[1]'] = k4[1]
    v['k4[2]'] = k4[2]
    return v

def legacy_steps(v=None):
    n = v['n']
    sigma = v['sigma']
    rho = v['rho']
    beta = v['beta']
    dt = v['dt']
    y = [v['y[0]'], v['y[1]'], v['y[2]']]
    k1 = [v['k1[0]'], v['k1[1]'], v['k1[2]']]
    k2 = [v['k2[0]'], v['k2[1]'], v['k2[2]']]
    k3 = [v['k3[0]'], v['k3[1]'], v['k3[2]']]
    k4 = [v['k4[0]'], v['k4[1]'], v['k4[2]']]
    _s = 1
    _e = n
    return max(0, int(_e) - int(_s) + 1)

def legacy_final(v):
    y = [v['y[0]'], v['y[1]'], v['y[2]']]
    y = [v['y[0]'], v['y[1]'], v['y[2]']]
    y = [v['y[0]'], v['y[1]'], v['y[2]']]
    return [y[0], y[1], y[2]]
