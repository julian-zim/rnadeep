import os
import numpy as np


def parse_alignment(ali_path, dbn_path, filename):
	with open(os.path.join(ali_path, filename + '.aln')) as file:
		ali = list()
		line = file.readline()
		while line != '':
			if line != 'CLUSTAL \n':
				ali.append(line.split()[1])
			line = file.readline()

	with open(os.path.join(dbn_path, filename + '.dbn')) as dbn_file:
		dbn = dbn_file.readlines()[1].split()[0]

	return ali, dbn


def parse_alignments(ali_directory, dbn_directory):
	alis = list()
	dbns = list()
	filenames = os.listdir(ali_directory)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists(os.path.join(dbn_directory, filename_base + '.dbn')):
			new_ali, new_dbn = parse_alignment(ali_directory, dbn_directory, filename_base)
			alis.append(new_ali)
			dbns.append(new_dbn)
		else:
			print('Warning: No dbn data found for \'' + filename + '\'.')

	return alis, dbns


def draw_ali_sets(ali_directory, dbn_directory, splits=None):
	if splits is None:
		splits = [1]
	assert sum(splits) <= 1

	alis, dbns = parse_alignments(ali_directory, dbn_directory)

	if not (len(alis) == len(dbns)):
		raise ValueError('Something about the input file is odd.')

	num = len(alis)

	a = np.arange(num)
	for s in splits:
		ids = np.random.choice(a, int(num * s), replace=False)
		yield [(alis[i], dbns[i]) for i in sorted(ids)]
		a = [i for i in a if i not in ids]
