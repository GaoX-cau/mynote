#!/usr/bin/env python3


import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np


def input_file2(filename, loc=0):
    f = open(filename, 'r')
    i = next(f)
    x, y = [], []        
    for i in f:
        i = i.split()
        pos = (int(i[1]) + int(i[2]) - 1) / 2
        x.append(pos-loc)
        y.append(float(i[4]))
    f.close()
    x_new = np.linspace(min(x), max(x), 50)
    f = interp1d(x, y, kind='quadratic')
    y_smooth=f(x_new)
    return x, y
    #return x_new, y_smooth

teosinte_b = input_file2('sw4a-b.100bp.windowed.pi')
maize_b = input_file2('G4a-B.100bp.windowed.pi')
teosinte_u = input_file2('sw4a-u.100bp.windowed.pi', 1698)
maize_u = input_file2('G4a-U.100bp.windowed.pi', 1698)

plt.plot(teosinte_u[0] + teosinte_b[0], teosinte_u[1] + teosinte_b[1], label='teosinte')
plt.plot(maize_u[0] + maize_b[0], maize_u[1] + maize_b[1], label='maize')
plt.xticks([-1000, 0, 1000, 2000, 3000], ['-1000', 'TSS', '1000', '2000', '3000'])
plt.xlabel('Position (bp)')
plt.ylabel('Nucleotide diversity (π)')
plt.legend()
plt.show()

#######

union_u = set(teosinte_u[0]) & set(maize_u[0])
teosinte_pi_u = []
for n, i in enumerate(teosinte_u[0]):
    if i in union_u:
        teosinte_pi_u.append(teosinte_u[1][n])
teosinte_pi_u = np.array(teosinte_pi_u)

maize_pi_u = []
for n, i in enumerate(maize_u[0]):
    if i in union_u:
        maize_pi_u.append(maize_u[1][n])
maize_pi_u = np.array(maize_pi_u)

union_b = set(teosinte_b[0]) & set(maize_b[0])
teosinte_pi_b = []
for n, i in enumerate(teosinte_b[0]):
    if i in union_b:
        teosinte_pi_b.append(teosinte_b[1][n])
teosinte_pi_b = np.array(teosinte_pi_b)

maize_pi_b = []
for n, i in enumerate(maize_b[0]):
    if i in union_b:
        maize_pi_b.append(maize_b[1][n])
maize_pi_b = np.array(maize_pi_b)

union_u = list(union_u)
union_u.sort()
for i in union_u:
    print(i, end='\t')
print()

pi2 = []
for i in maize_pi_u / teosinte_pi_u :
    pi2.append(i)

union_b = list(union_b)
union_b.sort()
for i in union_b:
    print(i, end='\t')
print()

for i in maize_pi_b / teosinte_pi_b:
    pi2.append(i)

union = union_u + union_b

plt.plot(union, pi2, 'r', label='maize / teosinte')
plt.ylabel('πm / πT')
plt.xticks([-1000, 0, 1000, 2000, 3000], ['-1000', 'TSS', '1000', '2000', '3000'])
plt.show()
