#!/usr/local/bin/python3
from local.src.gor import GOR, prof_parse
import sys 
import numpy as np
matrices_dir = sys.argv[1]
sequence = sys.argv[2]

profile = prof_parse(sequence)
model = GOR()
model.load(matrices_dir)
prediction = model.predict(profile, padding=True)
print(prediction)