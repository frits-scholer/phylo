#!/usr/bin/env python

p = [0.009,0.165,0.0008,0.205,0.396,0.450,0.641,0.781,0.9,0.993]
print p
import statsmodels.api as sm
    # multiple-testing correction using BH procedure
q = sm.stats.multipletests(p)[1]
print q


