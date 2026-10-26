# --- English ---
# Preliminary audit: refits the R2 Table 2 linear models on 00b's merge, with and without duplicated
# rows and accent mismatches. Prints coefficients only (E2, E5).
# Auditoría R2 (Claude, 2026-10-22). Ver notas/auditoria_codigo.md. Ejecutar: python 00c_auditoria_replica_tabla2.py
from pathlib import Path
import pandas as pd, numpy as np, os
M=pd.read_csv(str(Path(__file__).resolve().parents[1]/'output'/'00_auditoria'/'nodos_merge_R2.csv'))
def fit(sub,var):
    s=sub.dropna(subset=[var])
    X=np.column_stack([np.ones(len(s)),(s.group=='off'),(s.group=='on'),(s.hs_partition==1),(s.group=='off')&(s.hs_partition==1),(s.group=='on')&(s.hs_partition==1)]).astype(float)
    b,*_=np.linalg.lstsq(X,s[var].values,rcond=None); return np.round(b,3),len(s)
for net,test,var in [('undirected','Acciones','AoA'),('undirected','Acciones','logFreq'),('directed','Fpalabra','logFreq'),('undirected','Animales','AoA')]:
    sub=M[(M.network==net)&(M.test==test)]
    ded=sub.drop_duplicates(subset=['network','test','group','name'])
    print(net,test,var,'con duplicados',fit(sub,var),'| sin duplicados (1ª coincidencia)',fit(ded,var))
print('--- emulando que las palabras con tilde no casan (lectura UTF-8 de csv latin1 en R)')
na=M.name.astype(str).str.contains(r'[^\x00-\x7f]',na=False)
M2=M.copy(); M2.loc[na,['logFreq','motor','Imageability','Concreteness','Familiarity','AoA']]=np.nan
for net,test,var in [('undirected','Acciones','AoA'),('undirected','Acciones','logFreq'),('directed','Fpalabra','logFreq'),('undirected','Animales','AoA'),('directed','Fpalabra','Familiarity')]:
    sub=M2[(M2.network==net)&(M2.test==test)]
    print(net,test,var,fit(sub,var))
