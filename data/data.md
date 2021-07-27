# The input data

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
