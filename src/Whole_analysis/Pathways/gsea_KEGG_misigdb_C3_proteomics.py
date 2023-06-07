
import gseapy


gseapy.enrichr(gene_list="/home/flomik/Desktop/Code-PHD/COCOMO_txn/processing/proteins_C3_co_expression.txt",description='C3_msigdb',gene_sets="/home/flomik/Desktop/Code-PHD/COCOMO_txn/data/GSEA/h.all.v7.5.1.symbols.gmt", outdir='/home/flomik/Desktop/Code-PHD/COCOMO_txn/results/GSEA', cutoff=0.5,verbose=True, background  = 2924)


gseapy.enrichr(gene_list="/home/flomik/Desktop/Code-PHD/COCOMO_txn/processing/proteins_C3_co_expression.txt",description='C3_KEGG',gene_sets="/home/flomik/Desktop/Code-PHD/COCOMO_txn/data/GSEA/c2.cp.kegg.v7.5.1.symbols_1_3_5.gmt", outdir='/home/flomik/Desktop/Code-PHD/COCOMO_txn/results/GSEA', cutoff=0.5,verbose=True, background  = 2924)








