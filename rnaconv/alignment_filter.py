import os
import sys
import numpy as np


default_max_single_diff = 0.1
default_max_doublet_diff = 0.1


def generate_equilibrium_frequencies(rfam_path, outpath):
	try:
		os.makedirs(os.path.join(outpath, 'frequencies',))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(outpath, 'frequencies/single',))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(outpath, 'frequencies/doublet'))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(outpath, 'frequencies/differences'))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(outpath, 'frequencies/differences/single',))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(outpath, 'frequencies/differences/doublet'))
	except FileExistsError:
		pass
	alidir_path = os.path.join(outpath, 'alignments')
	neidir_path = os.path.join(rfam_path, 'seed_neighbourhoods/nei')
	freqdir_path = os.path.join(rfam_path, 'seed_frequencies')

	base_to_ids = {
		'A': np.array([1, 0, 0, 0]),
		'C': np.array([0, 1, 0, 0]),
		'G': np.array([0, 0, 1, 0]),
		'U': np.array([0, 0, 0, 1]),
		'M': np.array([0.5, 0.5, 0, 0]),
		'R': np.array([0.5, 0, 0.5, 0]),
		'W': np.array([0.5, 0, 0, 0.5]),
		'S': np.array([0, 0.5, 0.5, 0]),
		'Y': np.array([0, 0.5, 0, 0.5]),
		'K': np.array([0, 0, 0.5, 0.5]),
		'V': np.array([1. / 3., 1. / 3., 1. / 3., 0]),
		'H': np.array([1. / 3., 1. / 3., 0, 1. / 3.]),
		'D': np.array([1. / 3., 0, 1. / 3., 1. / 3.]),
		'B': np.array([0, 1. / 3., 1. / 3., 1. / 3.]),
		'N': np.array([0.25, 0.25, 0.25, 0.25]),
	}

	for alifilename in os.listdir(alidir_path):
		filename_base = alifilename.split('.')[0]

		alignments = list()
		with open(os.path.join(alidir_path, alifilename)) as file:
			ali_idx = -1
			line = file.readline()
			while line != '':
				content = line.split()

				if content[1].isnumeric():
					ali_idx += 1
					alignments.append(list())
				else:
					alignments[ali_idx].append(content[1])

				line = file.readline()

		neighbourhood = set()
		with open(os.path.join(neidir_path, filename_base + '.nei')) as neifile:
			line = neifile.readline()
			while line != '':
				split_line = line.split()
				if len(split_line) < 3:
					split_line.append('-1')
				pair = tuple(sorted((int(split_line[1].split('|')[0]), int(split_line[2].split(':')[0]))))
				neighbourhood.add(pair)
				line = neifile.readline()

		with open(os.path.join(freqdir_path, 'single', filename_base + '.freq')) as single_freq_file:
			seed_single_frequencies = np.array([float(e) for e in single_freq_file.readline()[:-1].split()])
		with open(os.path.join(freqdir_path, 'doublet', filename_base + '.freq')) as doublet_freq_file:
			seed_doublet_frequencies = np.array([float(e) for e in doublet_freq_file.readline()[:-1].split()])

		single_frequencies_set = list()
		doublet_frequencies_set = list()
		single_differences_set = list()
		doublet_differences_set = list()
		for alignment in alignments:
			single_frequencies = np.ones(4)  # ones because of pseudocounts
			doublet_frequencies = np.ones(16)  # ones because of pseudocounts
			for seq in alignment:
				for pair in neighbourhood:
					if pair[0] == -1:
						char = seq[pair[1]]
						if char != '-':
							single_frequencies += base_to_ids[char]
					else:
						char_i, char_j = seq[pair[0]].upper(), seq[pair[1]].upper()
						if char_i != '-' and char_j != '-':
							summand = np.outer(base_to_ids[char_i], base_to_ids[char_j]).flatten()
							doublet_frequencies += summand

			# normalization
			sum_single_freq = sum(single_frequencies)
			if sum_single_freq > 0:
				single_frequencies /= sum_single_freq
			sum_double_freq = sum(doublet_frequencies)
			if sum_double_freq > 0:
				doublet_frequencies /= sum_double_freq

			single_frequencies_set.append(single_frequencies)
			doublet_frequencies_set.append(doublet_frequencies)
			single_differences_set.append(np.abs(single_frequencies - seed_single_frequencies))
			doublet_differences_set.append(np.abs(doublet_frequencies - seed_doublet_frequencies))

		with open(os.path.join(outpath, 'frequencies/single', filename_base + '.freq'), 'w') as outfile:
			outfile.write('\n'.join([' '.join([str(e) for e in single_frequencies]) for single_frequencies in single_frequencies_set]) + '\n')
		with open(os.path.join(outpath, 'frequencies/doublet', filename_base + '.freq'), 'w') as outfile:
			outfile.write('\n'.join([' '.join([str(e) for e in doublet_frequencies]) for doublet_frequencies in doublet_frequencies_set]) + '\n')
		with open(os.path.join(outpath, 'frequencies/differences/single', filename_base + '.freq'), 'w') as outfile:
			outfile.write('\n'.join([' '.join([str(e) for e in single_differences]) for single_differences in single_differences_set]) + '\n')
		with open(os.path.join(outpath, 'frequencies/differences/doublet', filename_base + '.freq'), 'w') as outfile:
			outfile.write('\n'.join([' '.join([str(e) for e in doublet_differences]) for doublet_differences in doublet_differences_set]) + '\n')


def filter_alignments(rfam_path, outpath, max_single_diff, max_doublet_diff):
	generate_equilibrium_frequencies(rfam_path, outpath)

	if max_single_diff is None:
		global default_max_single_diff
		max_single_diff = default_max_single_diff
	if max_doublet_diff is None:
		global default_max_doublet_diff
		max_doublet_diff = default_max_doublet_diff

	blacklist = set()

	for single_freq_filename in os.listdir(os.path.join(outpath, 'frequencies/differences/single')):
		with open(os.path.join(outpath, 'frequencies/differences/single', single_freq_filename)) as single_freq_file:
			differences = [float(e) for e in ' '.join(line[:-1] for line in single_freq_file.readlines()).split()]
		if not all([e < max_single_diff for e in differences]):
			blacklist.add(single_freq_filename.split('.')[0])

	for doublet_freq_filename in os.listdir(os.path.join(outpath, 'frequencies/differences/single')):
		with open(os.path.join(outpath, 'frequencies/differences/doublet', doublet_freq_filename)) as doublet_freq_file:
			differences = [float(e) for e in ' '.join(line[:-1] for line in doublet_freq_file.readlines()).split()]
		if not all([e < max_doublet_diff for e in differences]):
			blacklist.add(doublet_freq_filename.split('.')[0])

	for file in blacklist:
		os.remove(os.path.join(outpath, 'alignments', file + '.ali'))
		#os.remove(os.path.join(outpath, 'frequencies/single', file + '.freq'))
		#os.remove(os.path.join(outpath, 'frequencies/doublet', file + '.freq'))
		#os.remove(os.path.join(outpath, 'frequencies/differences/single', file + '.freq'))
		#os.remove(os.path.join(outpath, 'frequencies/differences/doublet', file + '.freq'))
		print('Removed alignment \'' + file + '\' due to having too high of a variance for the nucleotide frequencies')


def main():
	if not 3 <= len(sys.argv) <= 5:
		print('Usage: ./alignment_filter.py <rfam path> <outpath> [max single freq diff] [max doublet freq diff]')
		return -1
	rfam_path = sys.argv[1]
	outpath = sys.argv[2]

	max_single_diff = None
	if len(sys.argv) >= 4:
		max_single_diff = float(sys.argv[4])
	max_doublet_diff = None
	if len(sys.argv) == 5:
		max_doublet_diff = float(sys.argv[5])

	print('========== FILTERING SISSI ALIGNMENTS ==========')
	filter_alignments(rfam_path, outpath, max_single_diff, max_doublet_diff)
	print('Done.')


if __name__ == '__main__':
	main()
