"""Geometric charge sharing approximation; not a weighting-field calculation.
Gaussian width is a preset inspired by R015's ~4 mm FWHM charge profile.
Integrate over nearest-channel (Voronoi) cells, not the copper width alone.
No edge renormalization: out-of-array charge remains lost.
"""
import numpy as np
from scipy.special import ndtr
def fractions(centres,position,sigma):
 c=np.asarray(centres,dtype=float)
 if len(c)<2 or not np.all(np.diff(c)>0) or sigma<0 or not np.isfinite(position+sigma):raise ValueError('Invalid charge-sharing geometry')
 edges=np.r_[c[0]-(c[1]-c[0])/2,(c[:-1]+c[1:])/2,c[-1]+(c[-1]-c[-2])/2]
 if sigma==0:
  w=np.zeros(len(c));i=np.searchsorted(edges,position,side='right')-1
  if 0<=i<len(c):w[i]=1
  return w
 return np.diff(ndtr((edges-position)/sigma))
