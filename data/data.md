# The input data

1. 在 Linux Shell 用以下代码进行下载
2. 在浏览器中输入 url 进行下载，并存储在 data 文件夹

```sh
# for: pbmc3k_tutorial
# dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics
file="pbmc3k_filtered_gene_bc_matrices.tar.gz"
url="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
cd data
wget $url
tar -zxvf $file
cd ../
```

```sh
# for: spatial_vignette
# reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol
file="allen_cortex.rds"
url="https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1"
cd data
wget $url
cd ../
```

```sh
# for: spatial_vignette
# mouse single-cell RNA-seq hippocampus dataset, produced in Saunders*, Macosko*, et al. 2018
file="mouse_hippocampus_reference.rds"
url="https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=0"
cd data
wget $url
cd ../
```
