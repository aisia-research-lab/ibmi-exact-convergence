from pathlib import Path
import urllib.request, tarfile, shutil, tempfile, ssl
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"; DATA.mkdir(exist_ok=True)
names=[f"bcsstk{i:02d}" for i in range(1,9)]
bases=[
 "https://sparse-files.engr.tamu.edu/MM/HB/{name}.tar.gz",
 "https://suitesparse-collection-website.herokuapp.com/MM/HB/{name}.tar.gz",
]
for name in names:
    dest=DATA/f"{name}.mtx"
    if dest.exists():
        print("EXISTS",dest); continue
    ok=False
    for tmpl in bases:
        url=tmpl.format(name=name)
        try:
            print("GET",url)
            with tempfile.TemporaryDirectory() as td:
                tg=Path(td)/f"{name}.tar.gz"
                urllib.request.urlretrieve(url,tg)
                with tarfile.open(tg,"r:gz") as tf: tf.extractall(td)
                hits=list(Path(td).rglob(f"{name}.mtx"))
                if not hits: raise RuntimeError("archive contained no expected .mtx")
                shutil.copy2(hits[0],dest)
            print("OK",dest); ok=True; break
        except Exception as e:
            print("FAILED",url,repr(e))
    if not ok: print("MISSING",name)
