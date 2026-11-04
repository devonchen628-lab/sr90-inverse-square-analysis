from pathlib import Path
import csv,json,os,subprocess,hashlib,concurrent.futures,time
P=Path(__file__).resolve().parent;ROOT=P.parent.parent;runtime=Path('/Users/chenqiyue/Documents/Codex/2026-10-12/https-chatgpt-com-share-6acd35f5-5168/.runtime/geant4');src=P/'src';data=ROOT/'outputs/sr90_research/results/last_three.csv'
with data.open() as f:rows=list(csv.DictReader(f))
ds=[float(z['d_m']) for z in rows];assert len(ds)==32 and abs(ds[0]-.01)<1e-9 and abs(ds[-1]-.85)<1e-9
binary_air=ROOT/'work/offset21_shift31/build/task16_sensitivity.app/Contents/MacOS/task16_sensitivity';binary_vac=P/'build/far15_vac';env=os.environ.copy();env.pop('TASK16_GEOMETRY_VIEW',None);env.pop('TASK16_VIS_MACRO',None)
for k,n in [('G4LEDATA','EMLOW'),('G4LEVELGAMMADATA','PhotonEvaporation'),('G4RADIOACTIVEDATA','RadioactiveDecay'),('G4ENSDFSTATEDATA','ENSDFSTATE')]:env[k]=str(next((runtime/'share/Geant4/data').glob(n+'*')))
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
plan={'distances_m':ds,'offset_mm':18,'transverse_x_y_mm':[0,0],'source_recess_mm':6,'detector_recess_mm':3,'emitter_diameter_mm_assumed':22,'source_exit_diameter_mm':9.75,'detector_active_mm':[10,10,.3],'detector_opening_mm':[10,10],'threshold_keV':10,'air_batches':8,'air_source_histories_per_batch':400000,'air_replay_per_stage':16000,'vacuum_batches':12,'vacuum_source_histories_per_batch':1000000,'air_binary_sha256':sha(binary_air),'vacuum_binary_sha256':sha(binary_vac),'vacuum_source_sha256':sha(src/'geant_vacuum.cc'),'spectrum_sha256':sha(src/'spectra_unique.csv'),'uncertainty':'independent batch bootstrap; no geometry systematics','note':'same source and detector solids air/vac except medium; air sampling uses importance weighting and replay'}
(P/'plan.json').write_text(json.dumps(plan,indent=2));print('Distances',len(ds),flush=True)
def run(item):
 medium,i=item;seed=(2021000000 if medium=='air' else 2031000000)+i*100003;f=P/medium/f'batch_{i:02}.csv';f.parent.mkdir(exist_ok=True);log=f.with_suffix('.log');meta=f.with_suffix('.json');dtext=','.join(map(str,ds))
 if medium=='air':cmd=[str(binary_air),str(src/'spectra_unique.csv'),str(f),dtext,'400000','16000',str(seed),'11','4.875','smc','10','10','18','0','0']
 else:cmd=[str(binary_vac),str(src/'spectra_unique.csv'),str(f),dtext,'1000000',str(seed),'1']
 if meta.exists():assert json.loads(meta.read_text())['command']==cmd;print('Existing',medium,i,flush=True);return
 start=time.time()
 with log.open('w') as h:subprocess.run(cmd,env=env,stdout=h,stderr=subprocess.STDOUT,check=True)
 with f.open() as h:z=list(csv.DictReader(h))
 assert len(z)==6*len(ds),(medium,i,len(z))
 for x in z:
  if medium=='air':assert abs(float(x['active_distance_m'])-float(x['distance_m'])-.018)<1e-9 and abs(float(x['aperture_diameter_mm'])-9.75)<1e-9 and float(x['height_mm'])==0 and float(x['lateral_x_mm'])==0
  else:assert abs(float(x['offset_m'])-.018)<1e-9
 logtext=log.read_text();assert not any(e in logtext for e in ['G4Exception','Overlap is detected','Fatal','ERROR'])
 meta.write_text(json.dumps({'command':cmd,'elapsed_s':time.time()-start,'seed':seed},indent=2));print(medium,'batch',i+1,'complete',round(time.time()-start),'s',flush=True)
items=[('air',i) for i in range(8)]+[('vacuum',i) for i in range(12)]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(run,items))
