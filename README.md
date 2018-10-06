Markov Affinity-based Graph Imputation of Cells (MAGIC)
-------------------------------------------------------
Van Dijk, David, et al. "Recovering Gene Interactions from Single-Cell Data Using Data Diffusion." Cell (2018).

https://www.ncbi.nlm.nih.gov/pubmed/29961576

MAGIC has been implemented in Python3 and Matlab. See below for installation and Usage.

Overview and how to use MAGIC effectively
-------------------------------------------------------

MAGIC is an unsupervised non-parametric algorithm to impute and de-noise biological single-cell RNA-seq data sets. MAGIC achieves this objective by sharing information across similar cells via data diffusion. Algorithmically, MAGIC performs the following operations in the sequential order:

	1). Constructs a nearest-neighbor graph
	2). Converts the distance matrix to an affinity matrix, via Gaussian kernel
	3). Converts the affinity matrix to a Markov (or Transition) matrix, via row normalization
	4). Raises the Markov matrix to a power t, to emulate data diffusion
	5). The original data is then right-multiplied by the powered Markov matrix 

MAGIC has a few important parameters that play important role the quality of results. These are listed below:

1). *k* = This defines the number of nearest neighbor used to construct the graph. We recommend this to be small 	               enough that only local neighborhood of each cell is considered but big enough that that graph remains                       connected. By default, this is set to k = 30.

2). _ka_ = This dictates the standard deviation to be used in the Gaussian kernel. To elaborate, the standard 	 	                deviation in the Gaussian kernel for a given cell is set to be the distance to it's ka-th nearest neighbor.                  By default, this is set to _ka_ = _k_/3 = 10.

3). _t_ = This defines the power to which the Markov matrix is to be raised. This is arguably the most important 		  parameter. A very high _t_ can lead to over-smoothed results while a low _t_ can lead to noisy results. 	               The choice of _t_ can be context-dependent. It is crucial to spend time looking at the data and                             selecting _t_ appropriately. Typically, we recommend _t_ between 3 and 8.

For example Python notebooks, please go inside the notebooks directory and open: Magic_single_cell_RNAseq.ipynb

For example MATLAB usage, please go inside the matlab directory and open: test_magic.m

For a tutorial on how to use the GUI, please go inside docs and open: magic_tutorial.pptx


#### Installation and dependencies for the Python version
1. The Python3 version of MAGIC can be installed using:

        $> git clone git://github.com/pkathail/magic.git
        $> cd magic
        $> sudo -H pip3 install .

2. MAGIC depends on a number of `python3` packages available on pypi and these dependencies are listed in `setup.py`
All the dependencies will be automatically installed using the above commands

3. After pulling updates to MAGIC from github, the package must be uninstalled and reinstalled:
		
		$> sudo -H pip3 uninstall magic
		$> sudo -H pip3 install .
		
#### Usage

#### Instructions for the Matlab version
1. run_magic.m -- MAGIC imputation function
2. test_magic.m -- Shows how to run MAGIC. Also included is a function for loading 10x format data (load_10x.m)

##### Interactive command line
A tutorial on MAGIC usage and results visualization for single cell RNA-seq data can be found in this notebook: http://nbviewer.jupyter.org/github/pkathail/magic/blob/develop/notebooks/Magic_single_cell_RNAseq.ipynb


##### GUI
A python GUI is now available for MAGIC. After following the installation steps listed below, the GUI can be invoked using

        $> magic_gui.py

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
