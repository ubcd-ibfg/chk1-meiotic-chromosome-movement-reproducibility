"""FIXED VERSION of 02_series1_classifiers.py: folds in which train or test lacks one class are skipped (original crashed with an empty test set).
Supp. Fig. 2A–D: classifiers with cross-validation grouped by cell and cell-level permutation tests.
Runtime ~1-2 h on 2 cores (set NPERM lower for a quick check).
Output columns are those read by 05, 08 and 10 (g, acc_grouped, null_mean, p_perm; a, b, acc). BH q values are added in 08."""
from common import *
import itertools, warnings; warnings.filterwarnings('ignore')
from joblib import Parallel, delayed
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression; from sklearn.tree import DecisionTreeClassifier; from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier; from sklearn.svm import SVC; from sklearn.ensemble import RandomForestClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis; from sklearn.neural_network import MLPClassifier
NPERM=int(os.environ.get('NPERM',200))
d=pd.read_csv(os.path.join(OUT,'frames_series1.csv')); F=['area','minor','circularity','cal.speed','ang.speed']; d=d.dropna(subset=F)
MODELS={'Logistic regression':lambda s:LogisticRegression(max_iter=1000),'Decision tree':lambda s:DecisionTreeClassifier(random_state=s),'Naive Bayes':lambda s:GaussianNB(),
 'k-NN (k=5)':lambda s:KNeighborsClassifier(5),'SVM (RBF)':lambda s:SVC(),'Random forest':lambda s:RandomForestClassifier(100,random_state=s),
 'LDA':lambda s:LinearDiscriminantAnalysis(),'QDA':lambda s:QuadraticDiscriminantAnalysis(),'MLP':lambda s:MLPClassifier((64,32),max_iter=300,random_state=s)}
def cv_acc(X,y,grp,seed,model='Random forest',scale=True):
    """5-fold CV with all frames of a cell in the same fold; train and test folds balanced by undersampling."""
    r=np.random.default_rng(seed); acc=[]
    for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=seed).split(X,y,grp):
        def bal(idx):
            i0=idx[y[idx]==0]; i1=idx[y[idx]==1]; n=min(len(i0),len(i1)); return np.r_[r.choice(i0,n,False),r.choice(i1,n,False)]
        if min(np.bincount(y[tr],minlength=2))==0 or min(np.bincount(y[te],minlength=2))==0: continue   # FIX: degenerate fold (one class absent)
        tr,te=bal(tr),bal(te); Xtr,Xte=X[tr],X[te]
        if scale: sc=MinMaxScaler().fit(Xtr); Xtr,Xte=sc.transform(Xtr),sc.transform(Xte)
        m=MODELS[model](seed).fit(Xtr,y[tr]); acc.append((m.predict(Xte)==y[te]).mean())
    return np.mean(acc) if acc else np.nan
def data(a,b):
    x=d[d.genotype.isin([a,b])]; return x[F].values,(x.genotype==b).astype(int).values,x.particle.values
def perm(X,y,grp,k):
    r=np.random.default_rng(1000+k); cells=np.unique(grp); lab=dict(zip(cells,r.permutation([y[grp==c][0] for c in cells])))
    return cv_acc(X,np.array([lab[c] for c in grp]),grp,k,scale=False)
G=ORDER[1:]; res=[]
for g in G:
    X,y,grp=data('control',g); obs=np.mean([cv_acc(X,y,grp,s,scale=False) for s in range(3)])
    null=Parallel(n_jobs=-1)(delayed(perm)(X,y,grp,k) for k in range(NPERM))
    res.append(dict(g=g,acc_grouped=obs,null_mean=np.mean(null),p_perm=(1+np.sum(np.array(null)>=obs))/(NPERM+1))); print(res[-1],flush=True)
R=pd.DataFrame(res); R.to_csv(os.path.join(OUT,'suppfig2b_rf_vs_control.csv'),index=False)
B=pd.DataFrame(Parallel(n_jobs=-1)(delayed(lambda g,m,s:(g,m,s,cv_acc(*data('control',g),s,m)))(g,m,s) for g in G for m in MODELS for s in range(3)),columns=['genotype','classifier','seed','accuracy'])
B.groupby(['genotype','classifier']).accuracy.mean().unstack().rename_axis('g').to_csv(os.path.join(OUT,'suppfig2a_benchmark.csv'))
gini=[]
for g in G:
    X,y,_=data('control',g); m=RandomForestClassifier(100,random_state=0).fit(MinMaxScaler().fit_transform(X),y); gini.append([g]+list(m.feature_importances_))
pd.DataFrame(gini,columns=['g']+F).sort_values('g').to_csv(os.path.join(OUT,'suppfig2c_gini.csv'),index=False)
pw=[(a,b,np.mean(Parallel(n_jobs=-1)(delayed(cv_acc)(*data(a,b),s,'Random forest',False) for s in range(4)))) for a,b in itertools.combinations(ORDER,2)]
pd.DataFrame(pw,columns=['a','b','acc']).to_csv(os.path.join(OUT,'suppfig3_pairwise.csv'),index=False)
print('done')
