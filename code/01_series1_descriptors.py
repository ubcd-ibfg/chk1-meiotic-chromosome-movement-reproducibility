"""Fig. 1c, 1e, 1f, 1g; Supplementary Fig. 1; Supplementary Table 1 (Series 1 rows).
Input : horsetail_data.csv (frame-level ChroMo output). Unit of analysis: the meiocyte."""
from common import *
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
h=pd.read_csv(FRAMES,low_memory=False)
for c in ['cal.speed','ang.speed','cum.dist']: h[c]=pd.to_numeric(h[c].astype(str).str.rstrip(';'),errors='coerce')
d=h[h.label2.isin(S1)].copy(); d['genotype']=d.label2.map(S1)
d.to_csv(os.path.join(OUT,'frames_series1.csv'),index=False)
D6=['area','major','minor','circularity','convexity','cal.speed']
cell=d.groupby(['genotype','particle']).agg(**{k:(k,'mean') for k in D6},horsetail_start=('frame','min'),horsetail_end=('frame','max'),frames=('frame','size')).reset_index()
cell['post_horsetail_min']=-cell.horsetail_end            # frames are minutes relative to meiosis I onset (t = 0)
cell['horsetail_min']=cell.horsetail_end-cell.horsetail_start+1
per={}
for pid,x in d.groupby('particle'):                         # dominant period: FFT of y position, 4-10 min band
    x=x.sort_values('frame'); y=x['y'].values-x['y'].mean()
    if len(y)<30: continue
    f=np.fft.rfftfreq(len(y),1.0); P=np.abs(np.fft.rfft(y))**2; m=(f>=1/10)&(f<=1/4); per[pid]=1/f[m][np.argmax(P[m])]
cell['period_min']=cell.particle.map(per)
C=cell[cell.genotype=='control']
VARS={'cal.speed':'Mean CoM velocity (µm/min)','post_horsetail_min':'Post-horsetail duration (min)','period_min':'Dominant period (min)','area':'Area (µm²)',
      'horsetail_min':'Horsetail duration (min)','circularity':'Circularity','convexity':'Convexity','major':'Major axis (µm)','minor':'Minor axis (µm)'}
rows=[]
for v in VARS:
    ps=[mannwhitneyu(cell[cell.genotype==g][v].dropna(),C[v].dropna()).pvalue for g in ORDER[1:]]
    for g,p,q in zip(ORDER[1:],ps,bh(ps)): rows.append(dict(variable=v,genotype=g,median=cell[cell.genotype==g][v].median(),control_median=C[v].median(),P=p,q_BH=q))
stats=pd.DataFrame(rows); stats.to_csv(os.path.join(OUT,'stats_series1_vs_control.csv'),index=False)
# ---- Supplementary Fig. 1a: meiocyte size (separate measurement, transmitted light; input/SuppFig1a_cell_size.xlsx)
SIZE_FILE=os.path.join(ROOT,'raw_data','S_pombe_imaging','SuppFig1a_cell_size.xlsx')
STRAIN={663:'control',674:'exo1Δ',681:'snf22Δ',690:'mek1Δ',692:'cds1Δ',703:'rad3Δ',707:'tel1Δ',713:'rev3Δ',714:'rec12Δ',724:'fml1Δ',727:'dsk1Δ',732:'rdh54Δ',4213:'chk1Δ',1366:'rec8Δ'}
size=None
if os.path.exists(SIZE_FILE):
    raw=pd.read_excel(SIZE_FILE,header=None); recs=[]
    for j in range(raw.shape[1]):                     # wide block: header row = strain code, one column per strain
        h=raw.iloc[0,j]
        try: code=int(float(h))
        except (TypeError,ValueError):
            code={'rec8':1366,'chk1':4213}.get(str(h).strip().lower())   # columns labelled by gene name
            if code is None: continue
        if code not in STRAIN: continue
        for v in pd.to_numeric(raw.iloc[1:,j],errors='coerce').dropna(): recs.append((STRAIN[code],code,float(v)))
    size=pd.DataFrame(recs,columns=['genotype','strain','cell_size'])
    size.to_csv(os.path.join(OUT,'suppfig1a_cell_size.csv'),index=False)
    Cs=size[size.genotype=='control'].cell_size; gs_=[g for g in ORDER[1:] if g in set(size.genotype)]
    ps=[mannwhitneyu(size[size.genotype==g].cell_size,Cs).pvalue for g in gs_]
    sz=pd.DataFrame(dict(variable='cell_size',genotype=gs_,median=[size[size.genotype==g].cell_size.median() for g in gs_],control_median=Cs.median(),P=ps,q_BH=bh(ps)))
    sz.to_csv(os.path.join(OUT,'stats_suppfig1a_cell_size.csv'),index=False)
Z=pd.DataFrame({k:[(cell[cell.genotype==g][k].mean()-C[k].mean())/C[k].std() for g in ORDER[1:]] for k in D6},index=ORDER[1:])
Z.to_csv(os.path.join(OUT,'fig1f_zscores.csv'))
X=StandardScaler().fit_transform(cell[D6]); pca=PCA().fit(X); S=pca.transform(X); cell['PC1'],cell['PC2']=S[:,0],S[:,1]
pd.DataFrame(pca.components_[:2].T,index=D6,columns=['PC1','PC2']).assign(var_explained_PC1=pca.explained_variance_ratio_[0],var_explained_PC2=pca.explained_variance_ratio_[1]).to_csv(os.path.join(OUT,'fig1g_pca_loadings.csv'))
cell.to_csv(os.path.join(OUT,'cells_series1.csv'),index=False)
# Supplementary Table 1, Series 1 rows
st=d.groupby('genotype').agg(N_cells=('particle','nunique'),N_frames=('frame','size'),area=('area','mean'),area_sd=('area','std'),minor=('minor','mean'),minor_sd=('minor','std'),
   circ=('circularity','mean'),circ_sd=('circularity','std'),vel=('cal.speed','mean'),vel_sd=('cal.speed','std')).loc[ORDER]
st.to_csv(os.path.join(OUT,'supp_table1_series1.csv'))
# ---- figures
def fp(q): return '<0.0001' if q<1e-4 else (f'{q:.4f}' if q<0.001 else f'{q:.3f}')
def cellplot(ax,v,title):
    s=stats[stats.variable==v].set_index('genotype'); rng=np.random.default_rng(1)
    for i,g in enumerate(ORDER):
        x=cell[cell.genotype==g][v].dropna().values
        ax.boxplot(x,positions=[i],widths=0.55,showfliers=False,medianprops=dict(color='k',lw=1),boxprops=dict(lw=0.6),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
        ax.scatter(i+rng.uniform(-0.18,0.18,len(x)),x,s=7,color=CAT[g],alpha=0.85,lw=0,zorder=3)
        if i: q=s.loc[g,'q_BH']; ax.text(i,1.02,fp(q),transform=ax.get_xaxis_transform(),rotation=55,fontsize=5.5,ha='left',color='k' if q<0.05 else '#999')
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER,rotation=50,ha='right',fontstyle='italic'); ax.set_ylabel(VARS[v]); ax.set_title(title,pad=22)
fig=plt.figure(figsize=(9.5,7.2)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.25],hspace=0.75,wspace=0.3)
cellplot(fig.add_subplot(gs[0,0]),'cal.speed','c · velocity during horsetail'); cellplot(fig.add_subplot(gs[0,1]),'post_horsetail_min','e · post-horsetail duration')
LB=['Area','Major axis','Minor axis','Circularity','Convexity','Velocity']; ax=fig.add_subplot(gs[1,0]); im=ax.imshow(Z.values,cmap='RdBu_r',vmin=-4,vmax=4,aspect='auto')
for i in range(13):
    for j in range(6): ax.text(j,i,f'{Z.values[i,j]:.2f}',ha='center',va='center',fontsize=5.5,color='w' if abs(Z.values[i,j])>2.5 else 'k')
ax.set_xticks(range(6)); ax.set_xticklabels(LB,rotation=40,ha='right'); ax.set_yticks(range(13)); ax.set_yticklabels(ORDER[1:],fontstyle='italic')
plt.colorbar(im,ax=ax,shrink=0.7,label='z-score vs control'); ax.set_title('f · phenotypic deviation from control')
ax=fig.add_subplot(gs[1,1])
for g in ORDER: mk=(cell.genotype==g).values; ax.scatter(S[mk,0],S[mk,1],s=6,color=CAT[g],alpha=0.4,lw=0)
for j,l in enumerate(LB): ax.arrow(0,0,pca.components_[0,j]*2.5,pca.components_[1,j]*2.5,color='#999',lw=0.5,head_width=0.08); ax.text(pca.components_[0,j]*2.9,pca.components_[1,j]*2.9,l,fontsize=5.5,color='#777',ha='center')
for g in ORDER:
    mk=(cell.genotype==g).values; c=S[mk,:2].mean(0); ax.scatter(*c,s=38,color=CAT[g],edgecolor='k',lw=0.7,zorder=4)
    if g in ('rec8Δ','mek1Δ','rdh54Δ'): ax.annotate(g,c,xytext=(5,-8),textcoords='offset points',fontsize=6,fontstyle='italic')
    if g=='control': ax.annotate('control',c,xytext=(-30,-12),textcoords='offset points',fontsize=6)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f} %)'); ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f} %)'); ax.set_title('g · PCA of individual cells')
save(fig,'Fig1_c_e_f_g')
def sizeplot(ax):
    rng=np.random.default_rng(1); st=sz.set_index('genotype')
    for i,g in enumerate(ORDER):
        x=size[size.genotype==g].cell_size.values
        if len(x)==0:
            ax.text(i,0.03,'n.d.',transform=ax.get_xaxis_transform(),ha='center',fontsize=5.5,color='#999'); continue
        ax.boxplot(x,positions=[i],widths=0.55,showfliers=False,medianprops=dict(color='k',lw=1),boxprops=dict(lw=0.6),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
        ax.scatter(i+rng.uniform(-0.18,0.18,len(x)),x,s=7,color=CAT[g],alpha=0.85,lw=0,zorder=3)
        if i: q=st.loc[g,'q_BH']; ax.text(i,1.02,fp(q),transform=ax.get_xaxis_transform(),rotation=55,fontsize=5.5,ha='left',color='k' if q<0.05 else '#999')
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER,rotation=50,ha='right',fontstyle='italic'); ax.set_ylabel('Meiocyte size (µm)'); ax.set_title('a · meiocyte size',pad=22)
fig,axs=plt.subplots(4,2,figsize=(9.2,10.5)); axs=axs.ravel()
if size is not None: sizeplot(axs[0])
else: axs[0].axis('off'); axs[0].text(0.02,0.6,'a · meiocyte size: raw_data/SuppFig1a_cell_size.xlsx not found',fontsize=8,color='#C00000',transform=axs[0].transAxes)
for ax,(v,t) in zip(axs[1:],[('period_min','b · oscillation period'),('area','c · chromosome-mass area'),('horsetail_min','d · horsetail duration'),('circularity','e · circularity'),('convexity','f · convexity'),('major','g · major axis'),('minor','h · minor axis')]): cellplot(ax,v,t)
fig.tight_layout(); save(fig,'SuppFig1')
print('PCA variance explained:',(pca.explained_variance_ratio_[:2]*100).round(1)); print(stats[stats.q_BH<0.05][['variable','genotype','q_BH']].to_string(index=False))
