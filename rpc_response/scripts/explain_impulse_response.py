from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parents[2]/'deliverables/Garfield_RPC_Part1_Physics_20261021/新增波形图与数据'
plt.rcParams.update({'font.size':14,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False})
t=np.linspace(-10,150,3201)
def h(t,tau):
 u=np.maximum(t,0)/tau
 return np.where(t>=0,u*np.exp(1-u),0)
fig,ax=plt.subplots(1,2,figsize=(13,4.8))
for tau,c in zip([5,10,25],['#0000cc','#d60000','#008000']):
 ax[0].plot(t,h(t,tau),c=c,label=f'Preset shaping time {tau} ns')
ax[0].set(title='Unipolar shaping function',xlabel='Time after impulse (ns)',ylabel=r'$h_{\tau_s}(t)$',xlim=(-5,150),ylim=(-.03,1.08));ax[0].legend(fontsize=10)
t2=np.linspace(80,300,4401)
for name,delay,c in [('XR',103.3751807926,'#d60000'),('XL',216.6248192074,'#0000cc')]:
 ax[1].plot(t2,h(t2-delay,5),c=c,label=f'{name}: ideal delay {delay:.2f} ns')
ax[1].set(title='Assumed effective chain response',xlabel='Time after local impulse (ns)',ylabel=r'$Z_k(t)/G_k=h_{5\,ns}(t-\delta_k)$',xlim=(80,300),ylim=(-.03,1.08));ax[1].legend(fontsize=10)
for a in ax:a.minorticks_on();a.tick_params(which='both',top=True,right=True)
fig.tight_layout();fig.savefig(p/'冲激响应与单极成形函数.png',dpi=180);fig.savefig(p/'冲激响应与单极成形函数.pdf')
print(p/'冲激响应与单极成形函数.png')
