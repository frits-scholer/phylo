#!/usr/bin/env python
if __name__ == '__main__':
    from scipy import stats
    d1 = [0.001224, 0.002912, 0.002912, 0.003157, 0.003157, 0.003157]
    d2 = [0.001224, 0.002912, 0.002912]
    print stats.ks_2samp(d1, d2)
    
