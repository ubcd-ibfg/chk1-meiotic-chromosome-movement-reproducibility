"""Fig. 5 (S. cerevisiae TEL4L-GFP telomere tracking): panels A-E, G, H and Source Data.
Run:  python code/13_cerevisiae_analysis.py  (from any directory)
Inputs: raw_data/S_cerevisiae/{control_int5_PASS_excl12_filtered.csv, mutant_int5_PASS_excl42_filtered.csv};
feature definitions reused from original_scripts/cerevisiae_ML_analysis.py (with DT corrected to 1 min). Unit of analysis: the cell (46 control, 137 rdh54Δ)."""
import os as _os; _os.chdir(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))  # run from the Reproducibility folder
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold, cross_val_score, permutation_test_score
from sklearn.pipeline import Pipeline; from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression; from sklearn.tree import DecisionTreeClassifier; from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier; from sklearn.svm import SVC; from sklearn.ensemble import RandomForestClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
SRC=os.environ.get('CEREV_DIR','raw_data/S_cerevisiae'); OUT=os.environ.get('CEREV_OUT','output/S_cerevisiae/panels'); NPERM=int(os.environ.get('NPERM',200))
os.makedirs(OUT,exist_ok=True); DATA_OUT=os.environ.get('CEREV_DATA_OUT','output/S_cerevisiae'); os.makedirs(DATA_OUT,exist_ok=True)
plt.rcParams.update({'font.family':['Arial','Liberation Sans'],'font.size':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':0.7})
CC,MC='#2C6FB0','#E2524F'; PX=0.0663; DT=1.0
# ---- reuse the original feature extraction (cerevisiae_ML_analysis.py), with paths redirected
code=open(os.environ.get('CEREV_SCRIPT','original_scripts/S_cerevisiae/cerevisiae_ML_analysis.py')).read()
code=code.replace("/sessions/upbeat-awesome-wozniak/mnt/uploads/",SRC.rstrip('/')+'/')
exec(code[:code.index('# 3. ML CLASSIFICATION')])          # defines ctrl, mut, THRESH_85 (px), feat_all
DT=1.0   # min per frame (1-min time-lapse; the original ML script sets DT=5, which only rescales the MSD-fit intercept)
F=[c for c in feat_all.columns if c not in ('uid','label','n_frames')]
NAMES={'mean_step':'Mean step','std_step':'Step std. dev.','max_step':'Max step','p25_step':'Step 25th pct','p75_step':'Step 75th pct','iqr_step':'Step IQR',
 'kurt_step':'Step kurtosis','skew_step':'Step skewness','ng_param':'Non-Gaussianity','vacf1':'Velocity ACF lag 1','vacf2':'Velocity ACF lag 2','vacf3':'Velocity ACF lag 3',
 'rev_rate':'Reversal rate','alpha':'Diffusion exp. α','msd_ratio':'MSD ratio','msd_lag1':'MSD (lag 1)','rg':'Radius of gyration','span':'Span','conf':'Confinement',
 'large85':'N large steps (>85th)','large90':'N large steps (>90th)','large95':'N large steps (>95th)','frac_l85':'Frac. large steps (>85th)','frac_l95':'Frac. large steps (>95th)',
 'return_2fr':'Returns ≤2 min','return_3fr':'Returns ≤3 min','return_4fr':'Returns ≤4 min','rev_large':'Reversals after large step','acf_step_lag1':'Step-size ACF lag 1','acf_step_lag3':'Step-size ACF lag 3','msd_lag5':'MSD (lag 5)'}
# ---- per-cell table in µm (panels A, D, E)
TH=THRESH_85*PX; rows=[]
for lab,df in (('Control',ctrl),('rdh54Δ',mut)):
    for uid,g in df.groupby('uid'):
        g=g.sort_values('frame'); x=g.x_rel_px.values*PX; y=g.y_rel_px.values*PX; st=np.hypot(np.diff(x),np.diff(y))
        rows.append(dict(cell=uid,condition=lab,n_frames=len(g),mean_velocity_um_min=st.mean()/DT,
             msd_lag5_um2=np.mean((x[5:]-x[:-5])**2+(y[5:]-y[:-5])**2),n_large_steps=int((st>TH).sum()),any_large_step=bool((st>TH).any())))
cells=pd.DataFrame(rows); fa=feat_all.set_index('uid'); cells['step_acf_lag3']=cells.cell.map(fa['acf_step_lag3'])
cells.to_csv(os.path.join(DATA_OUT,'source_data_fig5_ADE_per_cell.csv'),index=False)
C=cells[cells.condition=='Control']; M=cells[cells.condition!='Control']
def mw(v): return stats.mannwhitneyu(C[v],M[v]).pvalue
def fp(p): return '<0.0001' if p<1e-4 else (f'{p:.4f}' if p<0.001 else f'{p:.3f}')
def save(fig,n):
    for e in ('pdf','svg'): fig.savefig(os.path.join(OUT,f'{n}.{e}'),bbox_inches='tight')
    plt.close(fig)
def violin(ax,v,yl,title=None):
    rng=np.random.default_rng(1)
    for i,(d,c) in enumerate(((C[v],CC),(M[v],MC))):
        p=ax.violinplot(d,[i],widths=0.75,showextrema=False); p['bodies'][0].set_facecolor(c); p['bodies'][0].set_alpha(0.45); p['bodies'][0].set_edgecolor('none')
        ax.boxplot(d,positions=[i],widths=0.18,showfliers=False,medianprops=dict(color='k',lw=1),boxprops=dict(lw=0.6),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
        ax.scatter(i+rng.uniform(-0.12,0.12,len(d)),d,s=4,color=c,alpha=0.6,lw=0,zorder=3)
    p=mw(v); ax.text(1,1.0,fp(p),transform=ax.get_xaxis_transform(),rotation=55,fontsize=6,color='k' if p<0.05 else '#999',ha='left',va='bottom')
    ax.set_xticks([0,1]); ax.set_xticklabels(['Control','rdh54Δ'],rotation=45,ha='right'); ax.get_xticklabels()[1].set_fontstyle('italic'); ax.set_ylabel(yl)
    if title: ax.set_title(title,pad=14)
# ---- A
fig,axs=plt.subplots(1,2,figsize=(3.0,2.2)); violin(axs[0],'mean_velocity_um_min','Velocity (µm/min)'); violin(axs[1],'msd_lag5_um2','MSD at lag 5 (µm²)')
fig.suptitle('TEL4L-GFP dynamics',fontsize=7.5,style='italic'); fig.tight_layout(); save(fig,'Fig5A_velocity_MSD')
# ---- B: balanced accuracy, 8 classifiers, repeated stratified 5-fold, permutation P
X=feat_all[F].fillna(0).values; y=feat_all.label.values
CLF={'Logistic regression':Pipeline([('s',StandardScaler()),('c',LogisticRegression(max_iter=2000))]),'Decision tree':DecisionTreeClassifier(random_state=42),
 'Naive Bayes':Pipeline([('s',StandardScaler()),('c',GaussianNB())]),'k-nearest neighbours':Pipeline([('s',StandardScaler()),('c',KNeighborsClassifier(5))]),
 'SVM (RBF)':Pipeline([('s',StandardScaler()),('c',SVC())]),'Random forest':RandomForestClassifier(200,random_state=42),
 'LDA':Pipeline([('s',StandardScaler()),('c',LinearDiscriminantAnalysis())]),'QDA':Pipeline([('s',StandardScaler()),('c',QuadraticDiscriminantAnalysis(reg_param=0.1))])}
rcv=RepeatedStratifiedKFold(n_splits=5,n_repeats=5,random_state=42); skf=StratifiedKFold(5,shuffle=True,random_state=42); res=[]
for n,c in CLF.items():
    ba=cross_val_score(c,X,y,cv=rcv,scoring='balanced_accuracy'); acc=cross_val_score(c,X,y,cv=rcv,scoring='accuracy')
    _,null,p=permutation_test_score(c,X,y,cv=skf,scoring='balanced_accuracy',n_permutations=NPERM,random_state=0,n_jobs=-1)
    res.append(dict(classifier=n,balanced_accuracy=ba.mean(),balanced_accuracy_sd=ba.std(),accuracy=acc.mean(),null_mean=null.mean(),P_perm=p)); print(res[-1],flush=True)
B=pd.DataFrame(res); B['majority_class_baseline_accuracy']=max(y.mean(),1-y.mean()); B.to_csv(os.path.join(DATA_OUT,'source_data_fig5B_classifiers.csv'),index=False)
fig,ax=plt.subplots(figsize=(3.3,2.4)); Bs=B.iloc[::-1].reset_index(drop=True)
for i,r in Bs.iterrows():
    s=r.P_perm<0.05; ax.barh(i,r.balanced_accuracy,xerr=r.balanced_accuracy_sd,color=MC if s else '#CFCFCF',height=0.65,error_kw=dict(lw=0.6,capsize=1.5))
    ax.text(r.balanced_accuracy+r.balanced_accuracy_sd+0.01,i,f'{r.balanced_accuracy:.2f}  P={r.P_perm:.3f}',va='center',fontsize=5.5,color='k' if s else '#999')
ax.axvline(0.5,ls='--',lw=0.6,color='#777'); ax.set_yticks(range(len(Bs))); ax.set_yticklabels(Bs.classifier); ax.set_xlim(0.3,0.85)
ax.set_xlabel('Balanced accuracy (5-fold CV × 5)'); ax.set_title('Classifier: rdh54Δ vs control',pad=6); save(fig,'Fig5B_classifiers_balanced')
# ---- C: RF importance (class-balanced RF) + univariate MW, top 13
rf=RandomForestClassifier(500,random_state=42,class_weight='balanced').fit(X,y)
imp=pd.DataFrame(dict(feature=F,name=[NAMES.get(f,f) for f in F],rf_importance=rf.feature_importances_,P_mannwhitney=[stats.mannwhitneyu(feat_all.loc[feat_all.label==0,f].dropna(),feat_all.loc[feat_all.label==1,f].dropna()).pvalue for f in F]))
imp['P_bonferroni']=(imp.P_mannwhitney*len(F)).clip(upper=1); imp=imp.sort_values('rf_importance',ascending=False); imp.to_csv(os.path.join(DATA_OUT,'source_data_fig5C_feature_importance.csv'),index=False)
t=imp.head(13).iloc[::-1].reset_index(drop=True); fig,ax=plt.subplots(figsize=(3.4,2.6))
for i,r in t.iterrows():
    s=r.P_mannwhitney<0.05; ax.barh(i,r.rf_importance,color='#F28E2B' if s else '#CFCFCF',height=0.7)
    ax.text(r.rf_importance+0.001,i,f'P={r.P_mannwhitney:.3f}'+('*' if s else ''),va='center',fontsize=5.5,color='#C05A00' if s else '#888')
ax.set_yticks(range(len(t))); ax.set_yticklabels(t.name); ax.set_xlabel('RF feature importance'); ax.set_title('RF importance and univariate P (uncorrected)',pad=6)
ax.set_xlim(0,t.rf_importance.max()*1.35); save(fig,'Fig5C_feature_importance')
# ---- D
fig,ax=plt.subplots(figsize=(1.6,2.2)); violin(ax,'step_acf_lag3','Step-size ACF at 3-min lag (a.u.)','Step-size ACF\nat 3-min lag'); ax.axhline(0,ls=':',lw=0.5,color='#999'); save(fig,'Fig5D_ACF_lag3')
# ---- E: % cells with ≥1 large step (>85th percentile of all steps), Wilson 95% CI, Fisher
def wilson(k,n,z=1.96):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
tab=[[int(C.any_large_step.sum()),len(C)],[int(M.any_large_step.sum()),len(M)]]
pE=stats.fisher_exact([[tab[0][0],tab[0][1]-tab[0][0]],[tab[1][0],tab[1][1]-tab[1][0]]])[1]
E_=pd.DataFrame([dict(condition=c,cells_with_large_step=k,cells=n,percent=100*k/n,CI_low=100*wilson(k,n)[0],CI_high=100*wilson(k,n)[1]) for c,(k,n) in zip(['Control','rdh54Δ'],tab)])
E_['P_fisher']=pE; E_['threshold_um']=TH; E_.to_csv(os.path.join(DATA_OUT,'source_data_fig5E_large_displacement.csv'),index=False)
fig,ax=plt.subplots(figsize=(1.5,2.2))
for i,r in E_.iterrows(): ax.bar(i,r.percent,0.6,color=[CC,MC][i],yerr=[[r.percent-r.CI_low],[r.CI_high-r.percent]],error_kw=dict(lw=0.6,capsize=2))
ax.text(1,1.0,fp(pE),transform=ax.get_xaxis_transform(),rotation=55,fontsize=6,ha='left',va='bottom'); ax.set_ylim(0,105); ax.set_ylabel(f'% cells with ≥1 step > {TH:.2f} µm')
ax.set_xticks([0,1]); ax.set_xticklabels(['Control','rdh54Δ'],rotation=45,ha='right'); ax.get_xticklabels()[1].set_fontstyle('italic'); ax.set_title('Large-displacement\nepisodes',pad=14); save(fig,'Fig5E_large_displacement')
# ---- G, H: representative cells (control C45 in field ctrol_3_2; rdh54Δ C36 in MAX_rdh54_02_1), first 25 frames
reps=[('Control','ctrol_3_2','C45',ctrl,CC),('rdh54Δ','MAX_rdh54_02_1','C36',mut,MC)]; G=[]; Hs=[]
fig,axs=plt.subplots(1,2,figsize=(3.6,1.9),sharex=True,sharey=True); fig2,axs2=plt.subplots(1,2,figsize=(3.6,1.5),sharey=True)
for k,(lab,field,cid,df,col) in enumerate(reps):
    g=df[(df.cell_id==cid)&(df.field_id.str.contains(field))].sort_values('frame').head(25)
    x=g.x_rel_px.values*PX; y=g.y_rel_px.values*PX; x-=x.mean(); y-=y.mean(); st=np.hypot(np.diff(x),np.diff(y))
    for i,(a,b) in enumerate(zip(x,y)): G.append(dict(condition=lab,cell=cid,frame=i+1,time_min=i*DT,x_um=a,y_um=b))
    for i,s_ in enumerate(st): Hs.append(dict(condition=lab,cell=cid,transition=i+1,step_um=s_,large=bool(s_>TH)))
    ax=axs[k]; ax.plot(x,y,color='#BBB',lw=0.5,zorder=1); sc=ax.scatter(x,y,c=np.arange(len(x)),cmap='viridis',s=10,zorder=2)
    ax.scatter(x[0],y[0],s=25,marker='o',facecolor='none',edgecolor='g',lw=1,zorder=3); ax.scatter(x[-1],y[-1],s=25,marker='s',facecolor='none',edgecolor='#D4A017',lw=1,zorder=3)
    ax.set_title(lab,style='italic' if k else 'normal'); ax.set_xlabel('X (µm)'); ax.set_aspect('equal')
    a2=axs2[k]; a2.bar(np.arange(1,len(st)+1),st,color=[MC if s_>TH else '#CCCCCC' for s_ in st],width=0.8); a2.axhline(TH,ls='--',lw=0.6,color='k')
    a2.set_title(lab,style='italic' if k else 'normal'); a2.set_xlabel('Frame')
axs[0].set_ylabel('Y (µm)'); fig.colorbar(sc,ax=axs,shrink=0.8,label='Frame'); save(fig,'Fig5G_trajectories')
axs2[0].set_ylabel('Step (µm)'); axs2[1].text(0.98,TH,f'threshold {TH:.2f} µm',transform=axs2[1].get_yaxis_transform(),ha='right',va='bottom',fontsize=5.5); fig2.tight_layout(); save(fig2,'Fig5H_steps')
pd.DataFrame(G).to_csv(os.path.join(DATA_OUT,'source_data_fig5G_trajectories.csv'),index=False); pd.DataFrame(Hs).to_csv(os.path.join(DATA_OUT,'source_data_fig5H_steps.csv'),index=False)
print('A: velocity P',round(mw('mean_velocity_um_min'),4),'MSD P',round(mw('msd_lag5_um2'),4),'| D: P',round(mw('step_acf_lag3'),4),'| E:',E_[['condition','percent','CI_low','CI_high']].round(1).values.tolist(),'P',round(pE,4))
print('min Bonferroni-adjusted univariate P:',round(imp.P_bonferroni.min(),3))
