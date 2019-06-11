# 简单转录组表达量分析流程

## Step1

```bash
sh run_rna-seq.sh
```

## Step2

对stringtie计算的表达量丰度进行PCA和层次聚类，可视化结果

```bash
./fpkm_cluster.py $(ls *.gene_abund.tab)
```

![PCA](PCA.png)

![hierarchy](hierarchy.png)

