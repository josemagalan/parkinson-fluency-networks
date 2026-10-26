# --- English ---
# Preliminary audit: rebuilds the R2 node/dictionary merge from the 18 R2 partition CSVs and reports
# row inflation from duplicated keys. Writes analisis/output/00_auditoria/nodos_merge_R2.csv (run 00a first).
# Auditoría R2 (Claude, 2026-10-22). Ver notas/auditoria_codigo.md. Ejecutar: python 00b_auditoria_merge_nodos.py
from pathlib import Path
import pandas as pd, glob, os, re, warnings; warnings.filterwarnings('ignore')
os.chdir(Path(__file__).resolve().parents[3]/'parkinson')
pat=re.compile('Acciones|Animales|Fpalabra|control|off|on|directed|undirected')
fr=[]
for net in ['directed','undirected']:
    for f in sorted(glob.glob(f'hub and spoke model/{net}/csvs/*.csv')):
        d=pd.read_csv(f,encoding='latin1'); info=pat.findall(f)
        d['test']=[i for i in info if i in ('Acciones','Animales','Fpalabra')][0]; d['group']=[i for i in info if i in ('control','off','on')][0]; d['network']=net
        fr.append(d)
N=pd.concat(fr,ignore_index=True)
print(N.groupby(['network','test','group']).agg(n=('name','size'),core=('hs_partition',lambda x:(x==0).sum())).unstack('group'))
lex=pd.read_excel('data/Lexesp data.xlsx',sheet_name='Hoja1'); lex=lex.rename(columns={lex.columns[0]:'w'})
aov=pd.read_excel('data/AoA_verbos.xlsx',sheet_name='Hoja1')[['Word','AoA_average']].set_axis(['w','AoA'],axis=1)
aon=pd.read_excel('data/AoA_sustantivos.xlsx',sheet_name='Hoja1')[['word','averageAoA']].set_axis(['w','AoA'],axis=1)
fr_=pd.read_excel('data/Buscador frecuencia oral Subtlex-Esp.xlsx',sheet_name='Subtlex-Esp')[['Word','Log freq.']].set_axis(['w','logFreq'],axis=1)
mo=pd.read_excel('data/buscador motor content.xlsx',sheet_name='motor_content').set_axis(['w','motor'],axis=1); mo['w']=mo.w.str.lower()
aoa=pd.concat([aov,aon])
M=N.copy()
for d,cols in [(fr_,['logFreq']),(mo,['motor']),(lex,['Imageability']),(lex,['Concreteness']),(lex,['Familiarity']),(aoa,['AoA'])]:
    M=M.merge(d[['w']+cols],left_on='name',right_on='w',how='left').drop(columns='w')
print('nodos',len(N),'filas tras merges',len(M))
infl=(M.groupby(['network','test','group']).size()/N.groupby(['network','test','group']).size()).unstack().round(3); print(infl)
lexd=lex[lex.w.duplicated(keep=False)]; print('lexesp duplicados:'); print(lexd[['w','Imageability','Familiarity']])
print('nodos con mayúsculas o tabuladores:', N.name.str.contains(r'[A-ZÁÉÍÓÚ]|\t| $').sum())
M.to_csv(str(Path(__file__).resolve().parents[1]/'output'/'00_auditoria'/'nodos_merge_R2.csv'),index=False)
# GLM replicate with statsmodels? 
try:
    import statsmodels.formula.api as smf
    sub=M[(M.network=='undirected')&(M.test=='Acciones')].dropna(subset=['AoA'])
    r=smf.ols('AoA ~ C(group, Treatment("control"))*C(hs_partition)',data=sub).fit(); print(r.params.round(3)); print('n',int(r.nobs))
except Exception as e: print('no statsmodels',e)
