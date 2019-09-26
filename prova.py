#!/usr/local/bin/python3
from local.src.myClass import GOR
import sys 
import numpy as np
matrices_dir = sys.argv[1]
np.save()
model = GOR()
model.load(matrices_dir)
print(model.dictionary['H'])