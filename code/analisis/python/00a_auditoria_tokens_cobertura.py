# --- English ---
# Preliminary audit of the R2 code (not part of the pipeline): tokens per participant, duplicated
# dictionary keys and per-token dictionary coverage. Reads ../parkinson/data; writes analisis/output/00_auditoria/.
# Auditoría R2 (Claude, 2026-10-22). Ver notas/auditoria_codigo.md. Ejecutar: python 00a_auditoria_tokens_cobertura.py
from pathlib import Path
import pandas as pd, glob, os, re
os.chdir(Path(__file__).resolve().parents[3]/'parkinson')
X=pd.read_excel('data/FLUIDEZ PARKINSON.xlsx',sheet_name=None)
task={'acciones park':'acciones','animales_park':'animales','palab F_park':'F','acciones cont':'acciones','animales_cont':'animales','palab F_Cont':'F'}
rows=[]
for s,df in X.items():
    for _,r in df.iterrows():
        ws=[w for w in r.iloc[2:].tolist() if isinstance(w,str) and w.strip()]
        for i,w in enumerate(ws): rows.append(dict(task=task[s],participant=r['participante'],group=r['grupo'],pos=i+1,word=w))
T=pd.DataFrame(rows)
print('tokens',len(T)); 
n=T.groupby(['task','group','participant']).size().rename('n').reset_index()
print(n.groupby(['task','group'])['n'].agg(['count','mean','std','min','max','sum']).round(2))
print('palabras con mayúscula/espacios:', T.word.str.contains(r'[A-ZÁÉÍÓÚ]|^\s|\s$').sum(), T.word[T.word.str.contains(r'[A-ZÁÉÍÓÚ]|^\s|\s$')].unique()[:15])
lex=pd.read_excel('data/Lexesp data.xlsx',sheet_name='Hoja1'); lex=lex.rename(columns={lex.columns[0]:'w'})
aov=pd.read_excel('data/AoA_verbos.xlsx',sheet_name='Hoja1')[['Word','AoA_average']].set_axis(['w','AoA'],axis=1)
aon=pd.read_excel('data/AoA_sustantivos.xlsx',sheet_name='Hoja1')[['word','averageAoA']].set_axis(['w','AoA'],axis=1)
fr=pd.read_excel('data/Buscador frecuencia oral Subtlex-Esp.xlsx',sheet_name='Subtlex-Esp')[['Word','Log freq.']].set_axis(['w','logFreq'],axis=1)
mo=pd.read_excel('data/buscador motor content.xlsx',sheet_name='motor_content').set_axis(['w','motor'],axis=1); mo['w']=mo.w.str.lower()
aoa=pd.concat([aov,aon])
for name,d in [('lexesp',lex),('aoa',aoa),('freq',fr),('motor',mo)]:
    print(name,'filas',len(d),'claves duplicadas',d.w.duplicated().sum())
D=T.copy()
D=D.merge(fr,left_on='word',right_on='w',how='left').drop(columns='w')
D=D.merge(mo,left_on='word',right_on='w',how='left').drop(columns='w')
D=D.merge(lex[['w','Imageability','Concreteness','Familiarity']],left_on='word',right_on='w',how='left').drop(columns='w')
D=D.merge(aoa,left_on='word',right_on='w',how='left').drop(columns='w')
print('filas tras merge',len(D),'(tokens',len(T),')')
for c in ['Imageability','Concreteness','Familiarity']: D[c]=pd.to_numeric(D[c],errors='coerce')
cov=D.groupby('task')[['logFreq','motor','Imageability','Concreteness','Familiarity','AoA']].apply(lambda g: g.notna().mean().round(3))
print('Cobertura por token (proporción no-NA):'); print(cov)
os.makedirs(str(Path(__file__).resolve().parents[1]/'output'/'00_auditoria'),exist_ok=True)
n.to_csv(str(Path(__file__).resolve().parents[1]/'output'/'00_auditoria'/'n_palabras_por_participante.csv'),index=False)
cov.to_csv(str(Path(__file__).resolve().parents[1]/'output'/'00_auditoria'/'cobertura_diccionarios_token.csv'))
# paired on/off n
p=n[n.group!='control'].pivot_table(index=['task','participant'],columns='group',values='n').reset_index()
print(p.groupby('task')[['on','off']].mean().round(2)); print('diff on-off mean', (p.on-p.off).groupby(p.task).mean().round(2))
