import os
import numpy as np


def parse_family(ali_filepath, dbn_filepath):
	"""
	Reads an alignment CLUSTAL file and neighbourhood Dot Bracket String file and combines them into a pair to be used
	for training.

		Parameters:
			ali_filepath (str): Path to the alignment CLUSTAL file
			dbn_filepath (str): Path to the neighbourhood Dot Bracket String file
	"""

	with open(ali_filepath) as file:
		ali = list()
		line = file.readline()
		while line != '':
			if line != 'CLUSTAL \n':
				ali.append(line.split()[1])
			line = file.readline()

	with open(dbn_filepath) as dbn_file:
		dbn = dbn_file.readlines()[1].split()[0]

	return ali, dbn


def parse_families(ali_dirpath, dbn_dirpath):
	"""
	Combines pairs of the same name of alignment CLUSTAL files and neighbourhood Dot Bracket String files found in the
	respective directories to be used for training.

		Parameters:
			ali_dirpath (str): Path to the directory containing the alignment CLUSTAL files
			dbn_dirpath (str): Path to the directory containing the neighbourhood Dot Bracket String files
	"""
	alis = list()
	dbns = list()
	filenames = os.listdir(ali_dirpath)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists(os.path.join(dbn_dirpath, filename_base + '.dbn')):
			new_ali, new_dbn = parse_family(os.path.join(ali_dirpath, filename),
											os.path.join(dbn_dirpath, filename_base + '.dbn'))
			alis.append(new_ali)
			dbns.append(new_dbn)
		else:
			print('Warning: No dbn data found for \'' + filename + '\'.')

	return alis, dbns


def draw_ali_sets(ali_directory, dbn_directory, splits=None):
	if splits is None:
		splits = [1]
	assert sum(splits) <= 1

	alis, dbns = parse_families(ali_directory, dbn_directory)

	if not (len(alis) == len(dbns)):
		raise ValueError('Something about the input file is odd.')

	num = len(alis)

	a = np.arange(num)
	for s in splits:
		ids = np.random.choice(a, int(num * s), replace=False)
		yield [(alis[i], dbns[i]) for i in sorted(ids)]
		a = [i for i in a if i not in ids]
