#!/usr/local/bin/python3

import os
import sys
import argparse
import numpy as np 
import pandas as pd

import magic

class ArgumentParserError(Exception):
	pass

class NewArgumentParser(argparse.ArgumentParser):
	def error(self, message):
		print(message)
		sys.exit(0)


def parse_args(args):
	p = NewArgumentParser(description='run MAGIC')
	p.add_argument('filetype', choices=['csv', '10x', '10x_HDF5', 'mtx'],
					help='what is the file type of your original data?')

	a = p.add_argument_group('data loading parameters')
	a.add_argument('-d', '--data-file', metavar='D', required=True,
					help='File path of input data file.')
	a.add_argument('-a', '--aff_mat_input_data_file', metavar='A', 
					help='File path of affinity matrix input data file.')
	a.add_argument('-o', '--output-file', metavar='O', required=True,
				   help='File path of where to save the MAGIC imputed data (in csv format).')
	a.add_argument('-q', '--output-file2', metavar='Q', required=True,
				   help='File path of where to save the diffusion map data (in csv format).')
	a.add_argument('-s', '--output-file3', metavar='S', required=True,
				   help='File path of where to save the markov affinity matrix (in csv format).')
	a.add_argument('-g', '--genome', metavar='G',
					help='Genome must be specified when loading 10x_HDF5 data.')
	a.add_argument('--gene-name-file', metavar='GN', 
					help='Gene name file must be specified when loading mtx data.')
	a.add_argument('--use-ensemble-ids', default=False, action='store_true',
					help='Use ensemble IDs instead of gene names.')
	a.add_argument('--cell-axis', metavar='CA', default='rows',
					choices=['rows', 'columns'], 
					help='When loading a csv, specify whether cells are on rows or columns (Default = \'rows\').')
	a.add_argument('--skip-rows', default=0, type=int,
					help='When loading a csv, number of rows to skip after the header row (Default = 0).')
	a.add_argument('--skip-columns', default=0, type=int,
					help='When loading a csv, number of columns to skip after the header columns (Default = 0).')

	n = p.add_argument_group('normalization/filtering parameters')
	n.add_argument('-n', '--no-normalize', default=True, action='store_false',
					help='Do not perform library size normalization on the data')
	n.add_argument('-l', '--log-transform', metavar='L', default=None, type=float,
					help='Log-transform data with the specified pseudocount.')
	n.add_argument('--mols-per-cell-min', default=0, type=int,
					help='Minimum molecules/cell to use in filtering.')
	n.add_argument('--mols-per-cell-max', default=np.inf, type=int,
					help='Maximum molecules/cell to use in filtering.')

	m = p.add_argument_group('MAGIC parameters')
	m.add_argument('-p', '--pca-components', metavar='P', default=20, type=int,
				   help='Number of pca components to use when running MAGIC (Default = 20).')
	m.add_argument('--pca-non-random', default=True, action='store_false',
				    help='Do not used randomized solver in PCA computation.')
	m.add_argument('-t', metavar='T', default=2, type=int,
					help='t parameter for running MAGIC (Default = 2).')
	m.add_argument('-k', metavar='K', default=9, type=int,
					help='Number of nearest neighbors to use when running MAGIC (Default = 9).')
	m.add_argument('-ka', metavar='KA', default=3, type=int,
					help='knn-autotune parameter for running MAGIC (Default = 3).')
	m.add_argument('-e', '--epsilon', metavar='E', default=1, type=int,
					help='Epsilon parameter for running MAGIC (Default = 1).')
	m.add_argument('-r', '--rescale', metavar='R', default=0, type=int,
					help='Percentile to rescale data to after running MAGIC (Default = 0).')

	w = p.add_argument_group('Diffusion Map parameters')
	w.add_argument('-c', '--n_diffusion_components', metavar='C', default=10, type=int,
					help='Number of diffusion map components to calculate (Default = 10).')

	try:
		return p.parse_args(args)
	except ArgumentParserError:
		raise


def main(args: list = None):
	args = parse_args(args)
	print(args)
	try:
		if args.filetype == 'csv':
			scdata = magic.mg.SCData.from_csv(os.path.expanduser(args.data_file), 
											  data_type='sc-seq', normalize=False, 
											  cell_axis=0 if args.cell_axis=='rows' else 1,
											  rows_after_header_to_skip=args.skip_rows,
											  cols_after_header_to_skip=args.skip_columns)
		elif args.filetype == 'mtx':
			scdata = magic.mg.SCData.from_mtx(os.path.expanduser(args.data_file),
											  os.path.expanduser(args.gene_name_file), normalize=False)
		elif args.filetype == '10x':
			scdata = magic.mg.SCData.from_10x(os.path.expanduser(args.data_file), normalize=False, 
											  use_ensemble_id=args.use_ensemble_ids)

		elif args.filetype == '10x_HDF5':
			scdata = magic.mg.SCData.from_10x_HDF5(os.path.expanduser(args.data_file), args.genome, 
												   normalize=False, use_ensemble_id=args.use_ensemble_ids)

		scdata.filter_scseq_data(filter_cell_min=args.mols_per_cell_min, filter_cell_max=args.mols_per_cell_max)

		if args.no_normalize:
			scdata = scdata.normalize_scseq_data()

		if args.log_transform != None:
			scdata.log_transform_scseq_data(pseudocount=args.log_transform)

		if args.aff_mat_input_data_file != None:
			aff_mat_input_data = pd.read_csv(args.aff_mat_input_data_file, sep=',')
		
			aff_mat_input_data.set_index(list(aff_mat_input_data.columns[[0]]), inplace=True)
				
			scdata.run_magic(aff_mat_input = aff_mat_input_data, n_pca_components=args.pca_components, random_pca=args.pca_non_random, t=args.t,
						 k=args.k, ka=args.ka, epsilon=args.epsilon, rescale_percent=args.rescale, n_diffusion_components=args.n_diffusion_components, output_file3 = args.output_file3)
		else: 
			scdata.run_magic(n_pca_components=args.pca_components, random_pca=args.pca_non_random, t=args.t,
						 k=args.k, ka=args.ka, epsilon=args.epsilon, rescale_percent=args.rescale, n_diffusion_components=args.n_diffusion_components, output_file3 = args.output_file3)

		scdata.magic.to_csv(os.path.expanduser(args.output_file))
		
		scdata.extended_data.to_csv(os.path.expanduser(args.output_file2))
		
		#scdata.diffusion_eigenvectors.to_csv(os.path.expanduser(args.output_file3))
		
		
		
	except:
		raise

if __name__ == '__main__':
	main(sys.argv[1:])
