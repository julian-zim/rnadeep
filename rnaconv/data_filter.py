import sys
import os
import numpy as np
import RNA


def obtain_and_compare_equilibrium_frequencies(ali_dirpath, neigh_dirpath,
											   orig_single_freq_dirpath, orig_doublet_freq_dirpath,
											   outpath):
	"""
	Extracts the equilibrium frequencies for unpaired single nucleotides and nucleotide pairs from the generated
	alignments and forms the differences to the already extracted equilibrium frequencies of the original alignments.

		Parameters:
			ali_dirpath (str): Path to the directory of the generated alignment files in CLUSTAL format
			neigh_dirpath (str): Path to the directory containing the alignment consensus structure files in dot bracket
				notation format
			orig_single_freq_dirpath (str): Path to the directory of the extracted single frequency files of the
				original alignments
			orig_doublet_freq_dirpath (str): Path to the directory of the extracted doublet frequency files of the
				original alignments
			outpath (str): Path to the directory in which to save the extracted unpaired single and paired nucleotide
				equilibrium frequencies and frequency differences
	"""

	os.makedirs(os.path.join(outpath, 'single'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'doublet'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'differences/single'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'differences/doublet'), exist_ok=True)

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

	for ali_filename in os.listdir(ali_dirpath):
		filename_base = ali_filename.split('.')[0]
		if not os.path.exists(os.path.join(neigh_dirpath, filename_base + '.dbn')):
			print('Warning: Couldn\'t obtain equilibrium frequencies of \"' + filename_base +
				  '\" since there is no neighbourhood information.')
			continue

		# read ali & dbn
		alignment = list()
		with open(os.path.join(ali_dirpath, ali_filename)) as alifile:
			alifile.readline()
			line = alifile.readline()
			while line != '':
				alignment.append(line.split()[1])
				line = alifile.readline()

		with open(os.path.join(neigh_dirpath, filename_base + '.dbn')) as dbnfile:
			dbnfile.readline()
			dbn = dbnfile.readline().split()[0]

		# convert dbn
		pairs = []
		stack = []
		for i, char in enumerate(dbn):
			if char == '(':
				stack.append(i)
			elif char == ')':
				if len(stack) != 0:
					j = stack.pop()
					pairs.append((j, i))
				else:
					raise ValueError('Mismatched parentheses neighbourhood file \'' + filename_base + '\'.')
			else:
				pairs.append((-1, i))
		if len(stack) != 0:
			raise ValueError('Mismatched parentheses neighbourhood file \'' + filename_base + '\'.')
		pairs = set(sorted(pairs))

		# generate frequencies
		single_frequencies = np.ones(4)  # ones because of pseudocounts
		doublet_frequencies = np.ones(16)  # ones because of pseudocounts
		for pair in pairs:
			for seq in alignment:
				if pair[0] == -1:
					char = seq[pair[1]]
					if char != '-':
						single_frequencies += base_to_ids[char]
				else:
					char_i, char_j = seq[pair[0]].upper(), seq[pair[1]].upper()
					if char_i != '-' and char_j != '-':
						summand = np.outer(base_to_ids[char_i], base_to_ids[char_j]).flatten()
						doublet_frequencies += summand

		# normalize
		sum_single_freq = sum(single_frequencies)
		if sum_single_freq > 0:
			single_frequencies /= sum_single_freq
		sum_double_freq = sum(doublet_frequencies)
		if sum_double_freq > 0:
			doublet_frequencies /= sum_double_freq

		# get original frequencies
		filename_seed = filename_base.split('_')[0]
		with open(os.path.join(orig_single_freq_dirpath, filename_seed + '.sfreq')) as sing_freq_file:
			single_seed_frequencies = np.array([float(value) for value in sing_freq_file.readline()[:-1].split()])
		with open(os.path.join(orig_doublet_freq_dirpath, filename_seed + '.dfreq')) as doub_freq_file:
			doublet_seed_frequencies = np.array([float(value) for value in doub_freq_file.readline()[:-1].split()])

		# get frequency differences
		sissi_frequencies = (single_frequencies, doublet_frequencies)
		seed_frequencies = (single_seed_frequencies, doublet_seed_frequencies)
		frequency_diffs = (sissi_frequencies[0] - seed_frequencies[0],
						   sissi_frequencies[1] - seed_frequencies[1])

		# save
		with open(os.path.join(outpath, 'single', filename_base + '.sfreq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in sissi_frequencies[0]]) + '\n')
		with open(os.path.join(outpath, 'doublet', filename_base + '.dfreq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in sissi_frequencies[1]]) + '\n')
		with open(os.path.join(outpath, 'differences/single', filename_base + '.sfreq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in frequency_diffs[0]]) + '\n')
		with open(os.path.join(outpath, 'differences/doublet', filename_base + '.dfreq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in frequency_diffs[1]]) + '\n')


def filter_data(ali_dirpath, seq_dirpath, neigh_dbn_dirpath, neigh_ct_dirpath, max_dbrs_deviation=20):
	"""
	Filters the data generated by the data_generator: In each alignment, every sequence is removed that
	deviates by over <max_dbrs_deviation> from the consensus structure that was used by SISSI to generate it. This is
	achieved by using RNAfold to predict the secondary structure and the base pair distance to compare it to the desired
	consensus structure. If this results in all sequences being removed, the whole alignment with the corresponding
	consensus structure and sequence is removed.
	Note: If families were generated, the consensus structures used for the alignment generation were generated by
		  RNAfold and saved into the neigh_dirpath.
		  If only alignments were generated, the consensus structures used for the generation were provided by the user,
		  most likely from a converted Rfam database, but then copied to the neigh_dirpath anyway, for the sake of
		  integrity.
		  Therefore, in both cases, neigh_dirpath can be used to retrieve the desired consensus structures to compare
		  the sequences with.

		Parameters:
			ali_dirpath (str): Path to the directory of the generated alignments.
			seq_dirpath (str): Path to the directory of the generated or copied sequences.
			neigh_dbn_dirpath (str): Path to the directory of the generated or copied dbn files.
			neigh_ct_dirpath (str): Path to the directory of the generated ct files.
			max_dbrs_deviation (int, None, optional): maximum allowed base pair distance deviation from the consensus structure in
				percent. Default is 20.
	"""

	if max_dbrs_deviation is None:
		max_dbrs_deviation = 20

	for ali_filename in os.listdir(ali_dirpath):
		filename = ali_filename.split('.')[0]

		with open(os.path.join(ali_dirpath, ali_filename), 'r') as alifile:
			ali = alifile.readlines()
		with open(os.path.join(neigh_dbn_dirpath, filename + '.dbn'), 'r') as dbnfile:
			cons_dbrs = dbnfile.readlines()[-1].split()[0]

		blacklist = list()
		for i, seq in enumerate([''.join(line.split()[1:]) for line in ali[1:]]):
			dbrs, _ = RNA.fold(seq)

			bp_diff = RNA.bp_distance(dbrs, cons_dbrs)
			if bp_diff / len(dbrs) > max_dbrs_deviation / 100:  # TODO: max bp distance should be len/2?
				blacklist.append(i)
				print('Deleting sequence ' + str(i) + ' of alignment ' + filename + ' due to its secondary structure '
									'predicted by RNAfold deviating '
									'by over ' + str(max_dbrs_deviation) + '%.')
		for i in reversed(blacklist):
			del ali[i+1]

		if len(ali) == 1:
			filename_base = filename.split('_')[0]
			if os.path.exists(os.path.join(seq_dirpath, filename_base + '.seq')):
				os.remove(os.path.join(seq_dirpath, filename_base + '.seq'))
			elif os.path.exists(os.path.join(seq_dirpath, filename + '.seq')):
				os.remove(os.path.join(seq_dirpath, filename + '.seq'))
			os.remove(os.path.join(ali_dirpath, ali_filename))
			os.remove(os.path.join(neigh_dbn_dirpath, filename + '.dbn'))
			os.remove(os.path.join(neigh_ct_dirpath, filename + '.ct'))
		else:
			with open(os.path.join(ali_dirpath, ali_filename), 'w') as alifile:
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

	# TODO: filter alignments who show a frequence distance aboe a certain threshold
	# or, you know, just dont allow the difference values to be higher or equal to 0.1

	#diffs = (list(), list())
	#diffs[0].append(read_frequencies(os.path.join(sissi_path, 'frequencies/differences/single', filename + '.freq')))
	#diffs[1].append(read_frequencies(os.path.join(sissi_path, 'frequencies/differences/doublet', filename + '.freq')))
	#diffs_flattened = ([vector for diffset in diffs[0] for vector in diffset],
	#					[vector for diffset in diffs[1] for vector in diffset])
	pass


def main():
	if not 2 <= len(sys.argv) <= 3:
		print('Usage: ./data_filter.py <data_directory_path> [<max_ss_deviation>]')
		return -1
	data_dirpath = sys.argv[1]
	max_dbrs_deviation = None
	if len(sys.argv) == 4:
		max_dbrs_deviation = sys.argv[3]

	ali_dirpath = os.path.join(data_dirpath, 'alignments')
	seq_dirpath = os.path.join(data_dirpath, 'sequences')
	neigh_dbn_dirpath = os.path.join(data_dirpath, 'neighbourhoods/dbn')
	neigh_ct_dirpath = os.path.join(data_dirpath, 'neighbourhoods/ct')
	if (not os.path.exists(ali_dirpath)
			or not os.path.exists(seq_dirpath)
			or not os.path.exists(neigh_dbn_dirpath)
			or not os.path.exists(neigh_ct_dirpath)):
		print('Warning: The given path doesn\'t look like a path to data generated by data_generator.py. We will filter'
			  'what we find.\n'
			  'Expected folder structure:\n'
			  '- data\n\t- alignments\n\t- neighbourhoods\n\t\t- ct\n\t\t- dbn\n\t- sequences')

	print('========== FILTERING DATA ==========')
	# obtain_and_compare_equilibrium_frequencies(ali_dirpath, neigh_dbn_dirpath,
	# 										   '../../../rfam/full/seed_frequencies/single',
	# 										   '../../../rfam/full/seed_frequencies/doublet',
	# 										   os.path.join(data_path, 'frequencies'))
	filter_data(ali_dirpath, seq_dirpath, neigh_dbn_dirpath, neigh_ct_dirpath, max_dbrs_deviation)
	print('Done.')


if __name__ == '__main__':
	main()
