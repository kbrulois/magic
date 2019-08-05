Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------
MAGIC is an imputation algorthm developed by Van Dijk, David, et al. "Recovering Gene Interactions from Single-Cell Data Using Data Diffusion." Cell (2018).

https://www.ncbi.nlm.nih.gov/pubmed/29961576

This is forked and modified from dpeerlab/magic. Removed support for the matlab version and the python gui. Added an R version and support for imputation across multiple batches. See below for installation and Usage.

Overview and how to use MAGIC effectively
-------------------------------------------------------
MAGIC is an unsupervised non-parametric algorithm to impyte and de-noise biological single-cell RNA-seq data sets. MAGIC achieves this objective by sharing information across similar cells via data diffusion. Algorithmically, MAGIC performs the following operations in the sequential order:

	1). Constructs a nearest-neighbor graph
	2). Converts the distance matrix to an affinity matrix, via a Gaussian kernel
	3). Converts the affinity matrix to a Markov (or Transition) matrix, via row normalization
	4). Raises the Markov matrix to a power t, to emulate data diffusion
	5). The original data is then right-multiplied by the powered Markov matrix
	
MAGIC has a few important parameters that play important rolw to determine the quality of results. They are listed below:

1). _k_ : This defines the number of nearest neighbors used to construct the graph. We recomment this to be small enough that only local neighborhood of each cell is considered but big enough that the graph remains connected. By default, this is set to k  = 30.

2). _ka_ : This dictates the standard deviation to be used in the Gaussian kernel. To elaborate, the standard deviation in the Gaussian kernel for a given cell is set to be the distance to it's ka-th nearest neighbor. By default, this is set to _ka_ = _k_/3 = 10.

3). _t_ : This defines the power to which the Markov matrix is to be raised. This is arguably the most important parameter. A very high _t_ can lead to over-smoothed results while a low _t_ can lead to noisy results. We provide an automatic way to detect _t_ in the paper. For this, we compute the degree of change between the imputed data at time _t_ and time _t-1_ and stop after this value stabilizes. With the increase in _t_, the data goes through a rapidly changing imputation regime followed by a smoothing regime. In the imputation regime, the diffusion learns the manifold structure and removes noise. At larger values of _t_, diffusing further would smooth out real biology. The knee-point determines an optimal _t_. However, the choice of _t_ can be context-dependent. While the automatic method can act as a good guide to choose optimal _t_, the underlying mathematics may not always recapitulate the true biology. Therefore, we always recommend to spend time looking at the data and to vet the values of _t_ accordingly. 

#### Installation and dependencies for the Python version
1. This Python3 version of MAGIC can be installed using:

        $> https://github.com/kbrulois/magic.git
        $> cd magic
        $> sudo -H pip3 install .

2. MAGIC depends on a number of `python3` packages available on pypi and these dependencies are listed in `setup.py`
All the dependencies will be automatically installed using the above commands

3. After pulling updates to MAGIC from github, the package must be uninstalled and reinstalled:
		
		$> sudo -H pip3 uninstall magic
		$> sudo -H pip3 install .
		

##### Command line script
MAGIC can be run using the command line script `MAGIC.py` with the following parameters:

		$> MAGIC.py -h
		usage: MAGIC.py [-h] -d D -o O [-g G] [--gene-name-file GN]
        		        [--use-ensemble-ids] [--cell-axis CA] [--skip-rows SKIP_ROWS]
                		[--skip-columns SKIP_COLUMNS] [-n] [-l L]
                		[--mols-per-cell-min MOLS_PER_CELL_MIN]
                		[--mols-per-cell-max MOLS_PER_CELL_MAX] [-p P]
                		[--pca-non-random] [-t T] [-k K] [-ka KA] [-e E] [-r R]
                		{csv,10x,10x_HDF5,mtx}
		
		run MAGIC

		positional arguments:
		  {csv,10x,mtx}         what is the file type of your original data?

		optional arguments:
		  -h, --help            show this help message and exit

		data loading parameters:
		  -d D, --data-file D   File path of input data file.
            -a A, --aff_mat_input_data_file A 
                                 File path of affinity matrix input data file.
		  -o O, --output-file O
		                        File path of where to save the MAGIC imputed data (in
		                        csv format).
		  -g G, --genome G      Genome must be specified when loading 10x_HDF5 data.
		  --gene-name-file GN   Gene name file must be specified when loading mtx
		                        data.
		  --use-ensemble-ids    Use ensemble IDs instead of gene names.
		  --cell-axis CA        When loading a csv, specify whether cells are on rows
		                        or columns (Default = 'rows').
		  --skip-rows SKIP_ROWS
		                        When loading a csv, number of rows to skip after the
		                        header row (Default = 0).
		  --skip-columns SKIP_COLUMNS
		                        When loading a csv, number of columns to skip after
		                        the header columns (Default = 0).
		
		normalization/filtering parameters:
		  -n, --no-normalize    Do not perform library size normalization on the data
		  -l L, --log-transform L
		                        Log-transform data with the specified pseudocount.
		  --mols-per-cell-min MOLS_PER_CELL_MIN
		                        Minimum molecules/cell to use in filtering.
		  --mols-per-cell-max MOLS_PER_CELL_MAX
		                        Maximum molecules/cell to use in filtering.

		MAGIC parameters:
		  -p P, --pca-components P
		                        Number of pca components to use when running MAGIC
		                        (Default = 20).
		  --pca-non-random      Do not used randomized solver in PCA computation.
		  -t T			t parameter for running MAGIC (Default = 6).
		  -k K			Number of nearest neighbors to use when running MAGIC
                        		(Default = 30).
		  -ka KA		knn-autotune parameter for running MAGIC (Default =
                        		10).
		  -e E, --epsilon E	Epsilon parameter for running MAGIC (Default = 1).
		  -r R, --rescale R	Percentile to rescale data to after running MAGIC
                        		(Default = 99).
##### Installation and dependencies for the R version
The R version can be installed using:
                        
                        install.packages("devtools")
                        devtools::install_github("kbrulois/magic", ref = "magicBatch")

