import sys
import os
import numpy as np
import RNA


default_max_dbrs_deviation = 20  # in percent


def obtain_sissi_frequencies(path, rfam_path):
	"""
	Extracts the equilibrium frequencies for unpaired single nucleotides and nucleotide pairs from the generated
	alignment and forms the differences to the already extracted equilibrium frequencies of the original rfam
	alignments.

		Parameters:
			path (str): path directory of the generated families. Expects the following folder structure:
				- alignments
				- neighbourhoods
					- ct
					- dbn
			rfam_path (str): path to the directory of the convered rfam database. Only used to retrieve the alignment,
			consensus structure and original equilibrium frequencies, therefore expects the following folder structure:
				- seed_alignments
				- seed_frequencies
					- single
					- doublet
				- seed_neighbourhoods
					- ct
					- dbn
					- nei
					- wuss
	"""

	os.makedirs(os.path.join(path, 'frequencies'), exist_ok=True)
	os.makedirs(os.path.join(path, 'frequencies/single'), exist_ok=True)
	os.makedirs(os.path.join(path, 'frequencies/doublet'), exist_ok=True)
	os.makedirs(os.path.join(path, 'frequencies/differences'), exist_ok=True)
	os.makedirs(os.path.join(path, 'frequencies/differences/single'), exist_ok=True)
	os.makedirs(os.path.join(path, 'frequencies/differences/doublet'), exist_ok=True)

	sissi_alidir_path = os.path.join(path, 'alignments')
	freqdir_path = os.path.join(rfam_path, 'seed_frequencies')
	neidir_path = os.path.join(rfam_path, 'seed_neighbourhoods/nei')

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

	for ali_filename in os.listdir(sissi_alidir_path):
		filename = ali_filename.split('.')[0]
		file_seedname = filename.split('_')[0]

		# get neighbourhoods
		neighbourhood = set()
		with open(os.path.join(neidir_path, file_seedname + '.nei')) as neifile:
			line = neifile.readline()
			while line != '':
				split_line = line.split()
				if len(split_line) < 3:
					split_line.append('-1')
				pair = tuple(sorted((int(split_line[1].split('|')[0]), int(split_line[2].split(':')[0]))))
				neighbourhood.add(pair)
				line = neifile.readline()

		# get sissi alignments
		with open(os.path.join(sissi_alidir_path, ali_filename), 'r') as alifile:
			sissi_alignment = [line.split()[1] for line in alifile.readlines()[1:]]

		# generate sissi frequencies
		single_frequencies = np.ones(4)  # ones because of pseudocounts
		doublet_frequencies = np.ones(16)  # ones because of pseudocounts
		for seq in sissi_alignment:
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

		sum_single_freq = sum(single_frequencies)
		if sum_single_freq > 0:
			single_frequencies /= sum_single_freq
		sum_double_freq = sum(doublet_frequencies)
		if sum_double_freq > 0:
			doublet_frequencies /= sum_double_freq

		sissi_frequencies = (single_frequencies, doublet_frequencies)

		# get seed frequencies
		with open(os.path.join(freqdir_path, 'single', file_seedname + '.freq')) as sing_freq_file:
			single_seed_frequencies = np.array([float(value) for value in sing_freq_file.readline()[:-1].split()])
		with open(os.path.join(freqdir_path, 'doublet', file_seedname + '.freq')) as doub_freq_file:
			doublet_seed_frequencies = np.array([float(value) for value in doub_freq_file.readline()[:-1].split()])

		seed_frequencies = (single_seed_frequencies, doublet_seed_frequencies)

		# get frequency differences
		frequency_diffs = (sissi_frequencies[0] - seed_frequencies[0],
						   sissi_frequencies[1] - seed_frequencies[1])

		# save
		with open(os.path.join(path, 'frequencies/single', filename + '.freq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in sissi_frequencies[0]]) + '\n')
		with open(os.path.join(path, 'frequencies/doublet', filename + '.freq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in sissi_frequencies[1]]) + '\n')
		with open(os.path.join(path, 'frequencies/differences/single', filename + '.freq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in frequency_diffs[0]]) + '\n')
		with open(os.path.join(path, 'frequencies/differences/doublet', filename + '.freq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in frequency_diffs[1]]) + '\n')


def filter_alignments(path, rfam_path, max_dbrs_deviation):
	"""
	Filters the alignments generated by the alignment_generator: In each alignment, every sequence is removed that
	deviates by over <max_dbrs_deviation> from the given consensus structure. This is achieved using RNAfold to predict
	the structure and the base pair distance to compare it to the desired consensus structure. If this results in all
	sequences being removed, the whole alignment file with the corresponding consensus structure file copy is removed.

		Parameters:
			path (str): path directory of the generated families. Expects the following folder structure:
				- alignments
				- neighbourhoods
					- ct
					- dbn
			rfam_path (str): path to the directory of the convered rfam database. Only used to retrieve the consensus
			structure, therefore expects the following folder structure:
				- seed_neighbourhoods
					- ct
					- dbn
			max_dbrs_deviation (int): maximum allowed base pair distance deviation from the consensus structure in
				percent. Default is 20.
	"""

	for ali_filename in os.listdir(os.path.join(path, 'alignments')):
		filename = ali_filename.split('.')[0]
		with open(os.path.join(rfam_path, 'seed_neighbourhoods/dbn', filename.split('_')[0] + '.dbn'), 'r') as dbnfile:
			cons_dbrs = dbnfile.readlines()[-1].split()[0]

		with open(os.path.join(path, 'alignments', ali_filename), 'r') as alifile:
			ali = alifile.readlines()

		blacklist = list()
		for i, seq in enumerate([''.join(line.split()[1:]) for line in ali[1:]]):
			dbrs, _ = RNA.fold(seq)

			bp_diff = RNA.bp_distance(dbrs, cons_dbrs)
			if bp_diff / len(dbrs) > max_dbrs_deviation / 100:
				blacklist.append(i)
				print('Deleting sequence ' + str(i) + ' of alignment ' + filename + ' due to its secondary structure '
									'predicted by RNAfold deviating '
									'by over ' + str(max_dbrs_deviation) + '%.')
		for i in reversed(blacklist):
			del ali[i+1]

		if len(ali) == 1:
			os.remove(os.path.join(path, 'alignments', ali_filename))
			os.remove(os.path.join(path, 'neighbourhoods/dbn', filename + '.dbn'))
			os.remove(os.path.join(path, 'neighbourhoods/ct', filename + '.ct'))
		else:
			with open(os.path.join(path, 'alignments', ali_filename), 'w') as alifile:
				alifile.write(''.join(ali))

		'''command = 'RNAalifold --noPS < ' + str(os.path.join(sissi_path, 'alignments', ali_filename))
		result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
		if len(result.stdout) == 0:
			raise RuntimeError('RNAalifold failed to predict the secondary structure for \'' + filename + '\':\n' + result.stderr)
		mfe_dbrs = result.stdout.split('\n')[1].split()[0]

		dbrs_diff = RNA.bp_distance(dbrs, mfe_dbrs)
		if dbrs_diff / len(dbrs) > max_dbrs_deviation / 100:
			os.remove(os.path.join(sissi_path, 'alignments', ali_filename))
			print('Removed alignment \'' + filename + '\' due to the secondary structure predicted by RNAalifold deviating '
									'by over ' + str(max_dbrs_deviation) + '%.')'''

	# TODO: filter alignments who show a frequence distance above a certain threshold
	# maybe just dont allow the difference values to be higher or equal to 0.1

	#diffs = (list(), list())
	#diffs[0].append(read_frequencies(os.path.join(sissi_path, 'frequencies/differences/single', filename + '.freq')))
	#diffs[1].append(read_frequencies(os.path.join(sissi_path, 'frequencies/differences/doublet', filename + '.freq')))
	#diffs_flattened = ([vector for diffset in diffs[0] for vector in diffset],
	#					[vector for diffset in diffs[1] for vector in diffset])


def main():
	if not 3 <= len(sys.argv) <= 4:
		print('Usage: ./generator_filter.py <path> <rfam path> [<max_ss_deviation>]')
		return -1
	path = sys.argv[1]
	rfam_path = sys.argv[2]
	max_dbrs_deviation = None
	if len(sys.argv) == 4:
		max_dbrs_deviation = sys.argv[3]

	print('========== FILTERING SISSI ALIGNMENTS ==========')
	#obtain_sissi_frequencies(rfam_path, sissi_path)
	filter_alignments(path, rfam_path, default_max_dbrs_deviation if max_dbrs_deviation is None else max_dbrs_deviation)
	print('Done.')


if __name__ == '__main__':
	main()
