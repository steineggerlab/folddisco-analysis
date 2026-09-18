"""paired.py <classes.tsv> <ref_dir> <dir>... : mean Sens@1FP and paired delta vs ref (95% bootstrap CI). Run inside the dir holding the result dirs."""
import sys, os, subprocess, numpy as np, pandas as pd
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{B}/scripts")
from mcsa_class_metrics import answer_classes, ranked_ids
classes = answer_classes(sys.argv[1])
def sens(d):
    out = {}
    for q, ans in classes.items():
        p = os.path.join(d, f"{q}.tsv")
        if not os.path.exists(p): continue
        n = 0
        for t in ranked_ids(p):
            if t not in ans: break
            n += 1
        out[q] = n / len(ans)
    return pd.Series(out)
ref = sens(sys.argv[2]); rng = np.random.default_rng(0)
print(f"{'config':28s} {'n':>3s} {'mean':>7s} {'dmean':>8s} {'ci_lo':>8s} {'ci_hi':>8s} win/loss")
for d in sys.argv[2:]:
    s = sens(d); j = pd.concat([ref, s], axis=1, keys=["r", "s"]).dropna()
    dd = (j.s - j.r).to_numpy()
    bs = [rng.choice(dd, len(dd)).mean() for _ in range(2000)]
    print(f"{os.path.basename(d):28s} {len(j):3d} {j.s.mean():7.4f} {dd.mean():+8.4f} {np.percentile(bs,2.5):+8.4f} {np.percentile(bs,97.5):+8.4f} {(dd>0).sum()}/{(dd<0).sum()}")
