"""Shared settings. All paths are relative to this package folder: raw_data/ (inputs) and output/ (results)."""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))   # the Reproducibility folder
RAW=os.path.join(ROOT,'raw_data','S_pombe_imaging','submission_data_AFA_2026'); FRAMES=os.path.join(ROOT,'raw_data','S_pombe_imaging','horsetail_data.csv')
OUT=os.path.join(ROOT,'output','S_pombe_imaging'); os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({'font.family':['Arial','Liberation Sans','DejaVu Sans'],'font.size':7,'axes.spines.top':False,'axes.spines.right':False,
                     'axes.linewidth':0.7,'pdf.fonttype':42,'svg.fonttype':'none'})
S1={'WT_2':'control','Δrec12':'rec12Δ','Δrdh54':'rdh54Δ','Δexo1':'exo1Δ','Δrev3':'rev3Δ','Δfml1':'fml1Δ','Δtel1':'tel1Δ','Δdsk1':'dsk1Δ',
    'Δrad3':'rad3Δ','Δcds1':'cds1Δ','Δmek1':'mek1Δ','Δchk1':'chk1Δ','Δsnf22':'snf22Δ','Δrec8':'rec8Δ'}
ORDER=list(S1.values())
CAT={'control':'#7F7F7F','rec12Δ':'#E8743B','rdh54Δ':'#D6338A','exo1Δ':'#D6338A','rev3Δ':'#D6338A','fml1Δ':'#D6338A','tel1Δ':'#6A5ACD','dsk1Δ':'#6A5ACD',
     'rad3Δ':'#6A5ACD','cds1Δ':'#6A5ACD','mek1Δ':'#6A5ACD','chk1Δ':'#6A5ACD','snf22Δ':'#2A9D6F','rec8Δ':'#D9A21B'}
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); prev=1.0
    for i in range(n-1,-1,-1): prev=min(prev,p[o[i]]*n/(i+1)); q[o[i]]=prev
    return q
def save(fig,name):
    for ext in ('pdf','svg'): fig.savefig(os.path.join(OUT,f'{name}.{ext}'),bbox_inches='tight')
    plt.close(fig)
