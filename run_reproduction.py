import argparse,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parent
def run(args): subprocess.run([sys.executable,*args],cwd=R,check=True)
ap=argparse.ArgumentParser()
ap.add_argument("--mode",choices=["smoke","paper"],default="smoke")
a=ap.parse_args()
run(["-m","pytest","-q"])
run(["experiments/covariance_validation.py"]+(["--quick"] if a.mode=="smoke" else []))
run(["experiments/ordering_stress.py"]+(["--quick"] if a.mode=="smoke" else []))
if a.mode=="paper":
    print("Public SPD reproduction requires SuiteSparse files. Run: python scripts/download_suitesparse.py")
    print("Then: python experiments/public_spd_validation.py")
run(["scripts/make_figures.py"])
