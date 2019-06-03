#!/usr/bin/env python3
# for Illumina data

import sys
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')


def get_ppi(libs, n):
    ppi = []
    for comb in combinations(libs, n):
        ss = set()
        for i in comb:
            ss = ss | data_set[i]
        ppi.append(len(ss))
    ppi = np.array(ppi)
    return ppi.mean(), ppi.std()


data_set = {}
for filename in sys.argv[1:]:
    f = open(filename, 'r')
    data_set[filename] = set()
    for i in f:
        i = i.split()
        if int(i[0]) >= 3:
            ss = i[1] + '\t' + i[2]
            data_set[filename].add(ss)
    f.close()

libs = data_set.keys()

x, y, yerr = [], [], []
for i in range(1, len(libs)+1):
    mean, std = get_ppi(libs, i)
    print(mean, std)
    x.append(i)
    y.append(mean)
    yerr.append(std)

plt.errorbar(x, y, yerr=yerr, fmt='o-')
plt.xlabel('No. of screens')
plt.ylabel('No. of protein interactions')
plt.savefig('errorbar.png', dpi=300)

