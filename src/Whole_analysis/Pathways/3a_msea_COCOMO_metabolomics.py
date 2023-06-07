## run in processing/msea

import gseapy

enr = gseapy.enrichr(gene_list="/home/flomik/Desktop/Code-PHD/COCOMO_txn/processing/com_3_metabolites.txt",description="met_C3",gene_sets="/home/flomik/Desktop/Code-PHD/COCOMO_txn/processing/Metabolon_met_2l.gmt",outdir='/home/flomik/Desktop/Code-PHD/COCOMO_txn/results/MSEA',cutoff=0.5,verbose=True, background = 877)



