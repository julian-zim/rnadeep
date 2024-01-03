import os
import numpy as np


def parse_alignment_set(ali_path, dbn_path, filename):
	with open(ali_path + filename + '.ali') as file:
		alis = list()
		ali_idx = -1
		line = file.readline()
		while line != '':
			content = line.split()

			if content[1].isnumeric():
				ali_idx += 1
				alis.append(list())
			else:
				alis[ali_idx].append(content[1])

			line = file.readline()

	with open(dbn_path + filename + '.dbn') as dbn_file:
		dbn = dbn_file.readlines()[1].split()[0]
	dbns = [dbn for _ in alis]

	return alis, dbns


def parse_alignments(ali_directory, dbn_directory):
	alis = list()
	dbns = list()
	filenames = os.listdir(ali_directory)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists(dbn_directory + filename_base + '.dbn'):
			new_alis, new_dbns = parse_alignment_set(ali_directory, dbn_directory, filename_base)
			alis += new_alis
			dbns += new_dbns
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
