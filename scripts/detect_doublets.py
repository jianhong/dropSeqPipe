import scrublet as scr
import scipy.io
import numpy as np
import os

input_dir = snakemake.input.input_dir
counts_matrix = scipy.io.mmread(input_dir + '/expression.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + 'features.tsv', delimiter='\t', column=0))
with open(input_dir + 'barcodes.tsv') as barcodes_file:
	barcodes = barcodes_file.readlines()

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=snakemake.config['DOUBLET_DETECTION']['min_counts'],
                                                          min_cells=snakemake.config['DOUBLET_DETECTION']['min_cells'], 
                                                          min_gene_variability_pctl=snakemake.config['DOUBLET_DETECTION']['min_gene_variability_pctl'], 
                                                          n_prin_comps=snakemake.config['DOUBLET_DETECTION']['n_prin_comps'])
with open(snakemake.output[0],'w') as results:
	results.write('cell,score\n')
	for cell,score in zip(barcodes,doublet_scores):
		results.write('{},{}\n'.format(cell.rstrip('-1\n'),score))