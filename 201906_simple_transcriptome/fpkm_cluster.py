#!/usr/bin/env python3


import sys
import pandas as pd
import numpy
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')


data = []
for filename in sys.argv[1:]:
    df = pd.read_csv(filename, sep='\t', index_col ='Gene ID')
    df = df.rename(columns={'TPM': filename})
    data.append(df[filename])

data = pd.concat(data, axis=1, join='inner')
X = data.iloc[ : , :].values.T

from sklearn.decomposition import PCA
pca = PCA(2, svd_solver='randomized')
pca.fit(X)
X_pca = pca.transform(X)
for i in range(len(X_pca)):
    plt.scatter(X_pca[i, 0], X_pca[i, 1], s=15, alpha=0.5, label=data.columns[i], linewidths=0)

plt.title('Principal component analysis of TPM')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend()
plt.tight_layout()
plt.savefig('PCA.png', dpi=300)
plt.close()

from scipy.cluster import hierarchy
Z = hierarchy.linkage(X, method ='ward')
hierarchy.dendrogram(Z, labels = data.columns)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('hierarchy.png', dpi=300)

