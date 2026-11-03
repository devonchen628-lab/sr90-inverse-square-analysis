"""Regenerate current analysis outputs without altering raw or archived files."""
from pathlib import Path
import hashlib, json, os, subprocess, sys, tempfile
P=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'sr90-mpl'))
sys.path.insert(0,str(P/'src'))
from comparison import analyze
from plotting import make_figure

def main():
    manifest=json.loads((P/'docs/provenance.json').read_text())
    raw_before={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (P/'data/raw').glob('*.csv')}
    for item in manifest:
        if item.get('destination','').startswith('data/raw/'):
            assert raw_before[Path(item['destination']).name]==item['sha256']
    for d in ['tables','figures']:(P/'results'/d).mkdir(parents=True,exist_ok=True)
    fits,data,sensitivity=analyze(P/'data/raw/Task16_2.csv',P/'results/tables')
    make_figure(fits,data,sensitivity,P/'results/figures')
    subprocess.run([sys.executable,str(P/'src/deadtime.py')],check=True)
    assert raw_before=={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (P/'data/raw').glob('*.csv')}
    from PIL import Image
    with Image.open(P/'results/figures/final_three_panel.png') as im:
        assert abs(im.info['dpi'][0]-300)<.1
        assert abs(im.width/300*25.4-178)<.1 and abs(im.height/300*25.4-130)<.1
    check=P/'results/tables/verification.json'
    record=json.loads(check.read_text());record.update(complete_pipeline_passed=True,raw_CSV_hashes_unchanged=True,figure_export_dimensions_checked=True)
    check.write_text(json.dumps(record,indent=2)+'\n')
    print('PASS: raw hashes, component sums, shared response, geometry, Jacobians, multistart and sensitivity checks.')

if __name__=='__main__':main()
