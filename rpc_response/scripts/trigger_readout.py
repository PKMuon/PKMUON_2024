"""Offline 20%-of-own-peak leading-edge timing, user preset, not hardware CFD.
Multiple significant peaks are flagged rather than automatically paired.
"""
import numpy as np
from scipy.signal import find_peaks
def trigger(t,signal,fraction=.2):
 t=np.asarray(t,float);y=np.asarray(signal,float)
 if len(t)!=len(y) or len(t)<3 or not np.all(np.diff(t)>0) or not np.isfinite(t).all() or not np.isfinite(y).all() or not 0<fraction<1:raise ValueError('Invalid waveform')
 peak=float(y.max());out={'time_ns':float('nan'),'peak':peak,'threshold':fraction*peak,'valid':False,'status':'no_positive_signal'}
 if peak<=0:return out
 if y[0]>=out['threshold']:out['status']='left_truncated';return out
 imax=int(np.argmax(y));cross=np.flatnonzero((y[:-1]<out['threshold'])&(y[1:]>=out['threshold'])&(np.arange(len(y)-1)<imax))
 if not len(cross):out['status']='no_rising_crossing';return out
 k=int(cross[0]);out['time_ns']=float(t[k]+(out['threshold']-y[k])*(t[k+1]-t[k])/(y[k+1]-y[k]))
 peaks,_=find_peaks(y,height=out['threshold'],prominence=.05*peak)
 out['significant_peaks']=len(peaks)
 if y[-1]>=out['threshold']:out['status']='right_truncated'
 elif len(peaks)>1 or len(cross)>1:out['status']='ambiguous_multi_pulse'
 else:out.update(valid=True,status='ok')
 return out
def coordinate(t1,t81,centres,total_delay=320.,tap_delay=4.):
 c=np.asarray(centres,float);tap=(t1-t81+total_delay)/(2*tap_delay)
 if not np.isfinite(tap) or not 0<=tap<=len(c)-1:return float('nan')
 return float(np.interp(tap,np.arange(len(c)),c)-c[len(c)//2])
