"""Motif-rate statistics (Fig. 2e, 2j–l; Fig. 3j, k; Fig. 4d–g; Supp. Fig. 4) from the official motif counts
(*_motifs_and_windows_summary.csv). Each strain is compared with the control decomposed in the same batch
(exact conditional test for two Poisson rates, cells as exposure).
Observation windows (12-frame windows with t < −100 min, summed over cells; column "windows" of the summary files) are reported
for strain and control, with windows per cell and events per window (descriptive only; tests use cells as exposure)."""
from common import *
import glob
from scipy.stats import binomtest, beta
def test(k1,n1,k0,n0):
    k1,k0=int(k1),int(k0)
    if k1+k0==0: return np.nan,np.nan,np.nan,np.nan
    p=binomtest(k1,k1+k0,n1/(n1+n0)).pvalue; K=k1+k0
    lo=beta.ppf(0.025,k1,k0+1) if k1>0 else 0.0; hi=beta.ppf(0.975,k1+1,k0) if k0>0 else 1.0
    f=lambda q: q/(1-q)*n0/n1 if q<1 else np.inf
    return ((k1/n1)/(k0/n0) if k0 else np.inf),f(lo),f(hi),p
rows=[]
for f in sorted(glob.glob(os.path.join(RAW,'*','*_motifs_and_windows_summary.csv'))):
    s=pd.read_csv(f); ctrl=s[s.motifC_mutant.isna()].iloc[0]; nC=int(ctrl.N_cells); wC=float(ctrl.windows)   # observation windows of the same-batch control
    for _,r in s.dropna(subset=['motifC_mutant']).iterrows():
        for M in 'AC':
            rr,lo,hi,p=test(r[f'motif{M}_mutant'],r.N_cells,r[f'motif{M}_control'],nC)
            rows.append(dict(summary_file=os.path.basename(f),strain=r.strain,motif=M,events=r[f'motif{M}_mutant'],cells=r.N_cells,rate=r[f'motif{M}_mutant']/r.N_cells,
                 control_events_same_batch=r[f'motif{M}_control'],control_cells=nC,control_rate=r[f'motif{M}_control']/nC,rate_ratio=rr,CI_low=lo,CI_high=hi,P=p,
                 windows=r.windows,control_windows=wC,windows_per_cell=r.windows/r.N_cells,control_windows_per_cell=wC/nC,
                 events_per_window=r[f'motif{M}_mutant']/r.windows,control_events_per_window=r[f'motif{M}_control']/wC))
T=pd.DataFrame(rows); T.to_csv(os.path.join(OUT,'motif_rates_vs_same_batch_control.csv'),index=False)
def get(file,strain,M):
    x=T[(T.summary_file==file)&(T.strain==strain)&(T.motif==M)].iloc[0]; return x
extra=[]
for a,b,file in [('rdh54_chk1','rdh54','Figure_rdh54_double_mutants_motifs_and_windows_summary.csv'),('rdh54_rad3','rdh54','Figure_rdh54_double_mutants_motifs_and_windows_summary.csv'),
                 ('rdh54_cds1','rdh54','Figure_rdh54_double_mutants_motifs_and_windows_summary.csv'),('rdh54_rec12','rdh54','Figure_rdh54_double_mutants_motifs_and_windows_summary.csv'),
                 ('hrs1_t16A_rdh54','rdh54','Figure_Microtubules_Actin_microtubules_motifs_and_windows_summary.csv'),('klp6_rdh54','rdh54','Figure_Microtubules_Actin_microtubules_motifs_and_windows_summary.csv')]:
    x,y=get(file,a,'C'),get(file,b,'C'); rr,lo,hi,p=test(x.events,x.cells,y.events,y.cells); extra.append(dict(comparison=f'{a} vs {b}',rate_ratio=rr,CI_low=lo,CI_high=hi,P=p))
ser=[('Figure_1_motifs_and_windows_summary.csv','rdh54'),('Figure_rdh54_double_mutants_motifs_and_windows_summary.csv','rdh54'),
     ('Figure_Rodrigo_motifs_and_windows_summary.csv','rdh54_hyg'),('Figure_Microtubules_Actin_latA_motifs_and_windows_summary.csv','rdh54_nolatA')]
xs=[get(f,s,'C') for f,s in ser]; K1=sum(x.events for x in xs); N1=sum(x.cells for x in xs); K0=sum(x.control_events_same_batch for x in xs); N0=sum(x.control_cells for x in xs)
rr,lo,hi,p=test(K1,N1,K0,N0); extra.append(dict(comparison=f'rdh54 pooled 4 series ({int(K1)}/{N1} vs {int(K0)}/{N0})',rate_ratio=rr,CI_low=lo,CI_high=hi,P=p))
pd.DataFrame(extra).to_csv(os.path.join(OUT,'motif_key_comparisons.csv'),index=False); print(pd.DataFrame(extra).round(4).to_string(index=False))
