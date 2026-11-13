"""Supp. Fig. 2b grey ticks: Random Forest accuracy when frames of the same cell may fall in both training and test folds
(ordinary stratified 5-fold CV on frames, 3 repetitions) – shows the optimism removed by grouping folds by cell."""
from common import *
import warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
d=pd.read_csv(os.path.join(OUT,'frames_series1.csv')); F=['area','minor','circularity','cal.speed','ang.speed']; d=d.dropna(subset=F)
def acc(X,y,seed):
    r=np.random.default_rng(seed); out=[]
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=seed).split(X,y):
        def bal(idx):
            i0=idx[y[idx]==0]; i1=idx[y[idx]==1]; n=min(len(i0),len(i1)); return np.r_[r.choice(i0,n,False),r.choice(i1,n,False)]
        tr,te=bal(tr),bal(te); m=RandomForestClassifier(100,random_state=seed).fit(X[tr],y[tr]); out.append((m.predict(X[te])==y[te]).mean())
    return np.mean(out)
rows=[]
for g in ORDER[1:]:
    x=d[d.genotype.isin(['control',g])]; X=x[F].values; y=(x.genotype==g).astype(int).values
    rows.append(dict(genotype=g,acc_ungrouped=np.mean([acc(X,y,s) for s in range(3)]))); print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUT,'suppfig2b_ungrouped.csv'),index=False)
