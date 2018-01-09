#!/usr/bin/env python
if __name__ == '__main__':
    import numpy as np
    np.random.seed(12345678)  #fix random seed to get the same result
    n1 = 200  # size of first sample
    n2 = 300  # size of second sample
    from scipy import stats
    d1 = stats.norm.rvs(size=n1, loc=0., scale=1)
    d2 = stats.norm.rvs(size=n2, loc=0.5, scale=1.5)
    d3 = stats.norm.rvs(size=n2, loc=0.01, scale=1.0)
    d4 = stats.norm.rvs(size=n2, loc=0.0, scale=1.0)
    print d1
    print d2
    print d3
    print d4
    print stats.ks_2samp(d1, d2)
    print stats.ks_2samp(d1, d3)
    print stats.ks_2samp(d1, d4)
    
