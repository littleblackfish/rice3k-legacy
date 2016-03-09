#! /usr/bin/env python2

from glob import glob,iglob
from numpy import load

for fname in iglob ('*.npz') :
    print  load(fname)['s']
    
