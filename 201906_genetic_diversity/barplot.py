#!/usr/bin/env python3


import matplotlib.pyplot as plt


maize = [0.0036944510600706703, 0.0070364841977869986]
teosinte = [0.020995801653868853, 0.01932037408822506]

x1 = [1,4]
x2 = [2,5]

plt.bar(x1, teosinte, label='Teosinte')
plt.bar(x2, maize, label='Maize')
plt.legend()
plt.xticks([1.5, 4.5], ['Upstream', 'Gene body'])
plt.ylabel('Nucleotide diversity (π)')
plt.xlim(0,6)
plt.show()

teosinte = [3.77914, 33.9781, 160.172]
maize = [4.38199, 27.7822, 94.0956]

x1 = [1, 4, 7]
x2 = [2, 5, 8]

plt.bar(x1, teosinte, label='Teosinte')
plt.bar(x2, maize, label='Maize')
plt.legend()
plt.xticks([1.5, 4.5, 7.5], ['Ear', 'Leaf', 'Stem'])
plt.ylabel('FPKM')
plt.xlim(0,9)
plt.show()