import numpy as np
import pandas as pd
import pyabf

abfs = [
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch1/18524003.abf",
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch4/18524006.abf",
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch6/18524012.abf",
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch8/18524018.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch2/18529011.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch4/18529019.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch5/18529021.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch7/18529026.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch1/20n02004.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch2/20n02005.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch5/20n02014.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch3/20n02005.abf",
]

traces =[]
file = 1

for s in abfs:
    abf = pyabf.ABF(s)
    x = abf.sweepX
    y = abf.sweepY
    traces.append(pd.DataFrame({'file': file, 'time': x, 'current': y}))
    file += 1

df = pd.concat(traces, ignore_index=True, sort=False)
df.to_csv("data/traces_for_noise_analysis.csv")