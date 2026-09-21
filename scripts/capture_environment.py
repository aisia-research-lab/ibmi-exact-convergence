from pathlib import Path
import json, platform, sys, os
import numpy, scipy, pandas, matplotlib
try:
    import psutil
except Exception:
    psutil=None
out={
  "python":sys.version,
  "platform":platform.platform(),
  "machine":platform.machine(),
  "processor":platform.processor(),
  "numpy":numpy.__version__,
  "scipy":scipy.__version__,
  "pandas":pandas.__version__,
  "matplotlib":matplotlib.__version__,
  "logical_cpus": os.cpu_count(),
  "memory_gb": round(psutil.virtual_memory().total/2**30,2) if psutil else None,
}
p=Path(__file__).resolve().parents[1]/"results"/"environment.json"
p.parent.mkdir(exist_ok=True)
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2)); print("WROTE",p)
