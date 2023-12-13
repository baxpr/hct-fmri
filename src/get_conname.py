#!/usr/bin/env python3

import pandas
import sys

confile = sys.argv[1]
connum = sys.argv[2]

cons = pandas.read_csv(confile, dtype=str)
conname = cons.ConName[cons.ConNum==connum].iloc[0]

print(conname)
