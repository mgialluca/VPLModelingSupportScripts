def planck(T,lam):
    '''
    Input
    -----
    T : K
    lam : m

    Output
    ------
    W/m2 per um per sr
    '''
    # SI per um
    import numpy as np
    c1 = 2*6.6260755e-34*(2.99792458e8)**2
    c2 = 6.6260755e-34*2.99792458e8/1.380658e-23
    return 1e-6 * c1 * lam**(-5) / (np.exp(c2/(T*lam)) - 1)
