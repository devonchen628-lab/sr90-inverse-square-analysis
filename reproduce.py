"""Regenerate final tables and figures; callable from any working directory."""
from pathlib import Path
import subprocess,sys,os,tempfile,shutil
P=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'sr90-mpl'))
raw=P/'data/raw/Task16_2.csv'
with tempfile.TemporaryDirectory(prefix='sr90-analysis-') as td:
 td=Path(td)
 for script,dest in [('models.py',td/'models'),('comparison.py',td/'comparison')]:
  subprocess.run([sys.executable,str(P/'src'/script),str(raw),'--output-dir',str(dest)],check=True)
  for f in dest.iterdir():
   if f.name=='data_used.csv':continue
   target=P/'results'/('figures' if f.suffix in ('.png','.pdf') else 'tables')/f.name
   shutil.copy2(f,target)
 subprocess.run([sys.executable,str(P/'src/deadtime.py')],check=True)
 import pandas as pd,numpy as np,json
 got=pd.read_csv(P/'results/tables/M0_M3_comparison.csv')
 expected=np.array([2894.5347284403992,1052.8962278268339,1022.7902136755895,216.60539762144833])
 assert np.allclose(got.total_deviance,expected,rtol=1e-7,atol=1e-4)
 assert abs(got.iloc[3].mu_per_m-.44077770817129375)<1e-6
 ref=pd.read_csv(P/'docs/reference_offset_sensitivity.csv')
 actual=pd.read_csv(P/'results/tables/offset_sensitivity_M0_M3.csv')
 cols=['total_deviance','normalization','mu_per_m']
 assert np.allclose(actual[cols],ref[cols],rtol=2e-6,atol=1e-6,equal_nan=True)
 (P/'results/tables/verification.json').write_text(json.dumps({'passed':True,'total_deviances':got.total_deviance.tolist(),'mu_per_m':float(got.iloc[3].mu_per_m),'offset_sensitivity_matches_reference':True,'python':sys.version},indent=2)+'\n')
