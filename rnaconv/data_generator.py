import random
import os
import subprocess
import argparse
import RNA


def generate_sequence_structure_pair(length=85, min_paired_sites_percent=20):
	"""
	Repeatedly generates a random sequence and predicts its secondary structure using RNAfold, until the structure
	has at least min_paired_sites paired sites.

		Parameters:
			length (int, optional): Length of the random sequence
			min_paired_sites_percent (int, optional): Minimal required sites to be paired in percent
	"""

	seq = ''
	dbrs = ''
	min_paired_sites = int(length * min_paired_sites_percent / 100)
	paired_sites = -1
	while paired_sites < min_paired_sites:
		seq = ''.join(random.choice(['A', 'C', 'G', 'U']) for _ in range(length))
		dbrs, _ = RNA.fold(seq)
		paired_sites = 0
		for char in dbrs:
			if char == '(' or char == ')':
				paired_sites += 1
	return seq, dbrs


def db_to_ct(dbn, seq):
	"""
	Converts the consensus structures contained in the dot bracket notation input file into the connect table format.

		Parameters:
			dbn (str): Secondary structure in dot bracket notation
			seq (str): Sequence
	"""

	ptable = list(RNA.ptable(dbn))

	column1 = list(range(1, ptable[0] + 1))
	column2 = [c for c in seq]
	column3 = [idx - 1 for idx in column1[:-1]] + [0]
	column4 = [idx + 1 for idx in column1]
	column5 = ptable[1:]
	column6 = column1

	content = list()
	content.append('{: >5}'.format(ptable[0]) + ' ENERGY =     0.0    1')
	for i in column3:
		content.append('{: >5}'.format(column1[i]) +
					   ' ' + column2[i] +
					   '{: >8}'.format(column3[i]) +
					   '{: >5}'.format(column4[i]) +
					   '{: >5}'.format(column5[i]) +
					   '{: >5}'.format(column6[i]) + '')

	return '\n'.join(content)


def generate_alignments(sissi_filepath, number,
						tree_filepath, neigh_filepath, sfreq_filepath, dfreq_filepath, ali_filepath,
						outpath):
	"""
	Generates <number> alternative alignments for the given tree-, consensus-structure-, single- & doublet-frequencies-
	and, optionally for readding indels, alignment-file, using:
		- RNAinverse to generate an ancestral sequence for the provided consensus structure
		- SISSI simulate homologous sequence alignments (taking the generated ancestral sequence, provided tree,
		  provided consensus structure and provided equilibrium frequencies as input).
	Note:
	The provided consensus structure will also be copied into an additional file per generated alignment, in order to
	create pairs of samples and tags to be used for training (during the process, the dbn files are converted to ct
	files, which are also saved).
	The generated ancestral sequence will also be saved to maintain integrity.

		Parameters:
			sissi_filepath (str): Path to the compiled sissi099 file
			number (int): The number of alignments to generate
			tree_filepath (str): Path to a directory containing tree files in the newick string format ('.seed_tree')
			neigh_filepath (str): Path to a directory containing neighbourhood files in the dot-bracket notation format ('.dbn')
			sfreq_filepath (str): Path to a directory containing files storing a single frequency vector ('.sfreq')
			dfreq_filepath (str): Path to a directory containing files storing a doublet frequency vector ('.dfreq')
			ali_filepath (str, None): Path to a directory containing alignment files in the clustal format ('.aln')
			outpath (str): The Path to which to write the generated sequences, alignments & copied consensus structures
	"""

	# check files
	filename_base = os.path.basename(tree_filepath).split('.')[0]
	if (not os.path.exists(tree_filepath)
			or not os.path.exists(neigh_filepath)
			or not os.path.exists(sfreq_filepath)
			or not os.path.exists(dfreq_filepath)
			or (ali_filepath is not None and not os.path.exists(ali_filepath))):
		print('Warning: Skipping \'' + filename_base + '\' as it is missing a tree, neighbourhood, frequency or'
													   'alignment file.')
		return

	# get equilibrium frequencies
	with open(sfreq_filepath, 'r') as single_freq_file:
		single_frequencies = single_freq_file.readline()[:-1]
	with open(dfreq_filepath, 'r') as doublet_freq_file:
		doublet_frequencies = doublet_freq_file.readline()[:-1]

	# get consensus structure & convert to ct
	with open(neigh_filepath, 'r') as dbn_file:
		# note: dbn files also save a sequence. If you used the rfam converter, this sequence is the consensus sequence
		# in the Rfam.seed stockholm file. It is not used for any generation, only passed through the file conversions.
		seq = dbn_file.readline()
		dbn = dbn_file.readline().split()[0]
	ct = db_to_ct(dbn, seq)

	# generate & save ancestral sequence
	#ancestral_sequence = ''.join(random.choice(['A', 'C', 'G', 'U']) for _ in range(len(dbn)))
	ancestral_sequence, _ = RNA.inverse_fold(None, dbn)

	seq_out_filepath = os.path.join(outpath, 'sequences', filename_base + '.seq')
	with open(seq_out_filepath, 'w') as file:
		file.write(ancestral_sequence + '\n')

	# save copies of the generated consensus structure (as dbn and ct)
	for i in range(number):
		# note: order of indexing not important since these are only copies of the same data
		filename = filename_base + '_a' + str(i)
		with open(os.path.join(outpath, 'neighbourhoods/dbn', filename + '.dbn'), 'w') as file:
			file.write(seq + dbn + ' (0.0)\n')
		with open(os.path.join(outpath, 'neighbourhoods/ct', filename + '.ct'), 'w') as file:
			file.write(ct + '\n')

	# generate alignments
	command = (sissi_filepath +
			   ' -k ' + str(seq_out_filepath) +
			   ' -nr' + str(neigh_filepath) +
			   ' -fs ' + single_frequencies +
			   ' -fd ' + doublet_frequencies +
			   ' -l' + str(len(dbn)) +
			   ' -a' + str(number) +
			   ' -oc' +
			   ' ' + str(tree_filepath))
	result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if len(result.stdout) == 0:
		print('Warning: Couldn\'t generate alignments for \'' + filename_base + '\':\n' + result.stderr)
		return

	# split generated alignment set into seperate alignments
	result_lines = result.stdout.split('\n')[:-1]  # for some reason a line without content is added in the end
	ali_sets = list()
	ali_idx = -1
	for line in result_lines:
		if line == 'CLUSTAL ':
			ali_sets.append(list())
			ali_idx += 1
		ali_sets[ali_idx].append(line)

	# re-add indels if wanted
	if ali_filepath is not None:
		with open(ali_filepath, 'r') as alifile:
			alifile.readline()
			line_idx = 1
			line = alifile.readline()
			while line != '':
				for i, char in enumerate(line.split()[1]):
					if char == '-':
						for ali in ali_sets:
							seq = list(ali[line_idx])
							seq[i + len(ali[line_idx].split()[0]) + 1] = '-'
							ali[line_idx] = ''.join(seq)

				line = alifile.readline()
				line_idx += 1

	# save alignments
	for idx, ali in enumerate(ali_sets):
		with open(os.path.join(outpath, 'alignments', filename_base + '_a' + str(idx) + '.aln'), 'w') as outfile:
			for line in ali:
				outfile.write(line + '\n')

	print('Successfully generated alignment' + ('s' if number > 1 else '') + ' for \'' + filename_base + '\'.')


def generate_alignment_set(sissi_filepath, number, tree_dirpath, neigh_dirpath, sfreq_dirpath, dfreq_dirpath,
						   ali_dirpath, outpath):
	"""
	Generates <number> alternative alignments for each tree file in the given tree-directory, searching in the
	respectively given consensus-structure-, single- & doublet-frequencies- and, optionally for readding
	ndels, alignment-directories for files of the same name to use.

		Parameters:
			sissi_filepath (str): Path to the compiled sissi099 file
			number (int): The number of alignments to generate
			tree_dirpath (str): Path to a directory containing tree files in the newick string format ('.seed_tree')
			neigh_dirpath (str): Path to a directory containing neighbourhood files in the dot-bracket notation format ('.dbn')
			sfreq_dirpath (str): Path to a directory containing files storing a single frequency vector ('.sfreq')
			dfreq_dirpath (str): Path to a directory containing files storing a doublet frequency vector ('.dfreq')
			ali_dirpath (str, None): Path to a directory containing alignment files in the clustal format ('.aln')
			outpath (str): The Path to which to write the generated sequences, alignments & copied consensus structures
	"""

	# check paths
	if not os.path.exists(sissi_filepath):
		raise FileNotFoundError('No sissi099 file found at \'' + sissi_filepath + '\'')
	for path in [tree_dirpath, neigh_dirpath, sfreq_dirpath, dfreq_dirpath]:
		if not os.path.exists(path):
			raise FileNotFoundError('No path named \'' + path + '\'')
	if ali_dirpath is not None and not os.path.exists(ali_dirpath):
		raise FileNotFoundError('No path named \'' + ali_dirpath + '\'')

	# make necessary directories
	os.makedirs(os.path.join(outpath, 'sequences'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/dbn'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/ct'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'alignments'), exist_ok=True)

	# generate
	filenames = os.listdir(tree_dirpath)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		generate_alignments(sissi_filepath,
							number,
							os.path.join(tree_dirpath, filename_base + '.seed_tree'),
							os.path.join(neigh_dirpath, filename_base + '.dbn'),
							os.path.join(sfreq_dirpath, filename_base + '.sfreq'),
							os.path.join(dfreq_dirpath, filename_base + '.dfreq'),
							os.path.join(ali_dirpath, filename_base + '.aln') if ali_dirpath is not None else None,
							os.path.join(outpath))
	print('Done.')


def generate_families(sissi_filepath, number, min_length, max_length,
					  tree_filepath, sfreq_filepath, dfreq_filepath,
					  outpath):
	"""
	Generates <number> families for the given tree- and single- & doublet-frequencies-file, using:
		- random ancestral sequences of uniformly distributed lengths up to <maxlength>
		- RNAfold to predict secondary structures for these sequences to be used as consensus structures for the
		  alignment generation
		- SISSI simulate corresponding homologous sequence alignments (taking the random ancestral sequences, provided
		  tree, predicted consensus structures, and provided equilibrium frequencies as input).
	Note:
	During the process, the generated dbn files are converted to ct files, which are also saved.

		Parameters:
			sissi_filepath (str): Path to the compiled sissi099 file
			number (int): The number of families to generate
			min_length (int): Minimum allowed length of the ancestral sequences used to generate the families
			max_length (int): Maximum allowed length of the ancestral sequences used to generate the families
			tree_filepath (str): Path to a tree file in the newick string format ('.seed_tree')
			sfreq_filepath (str): Path to a file containing a single frequency vector ('.sfreq')
			dfreq_filepath (str): Path to a file containing a doublet frequency vector ('.dfreq')
			outpath (str): The path to which to write the generated families
	"""

	# check files
	filename_base = os.path.basename(tree_filepath).split('.')[0]
	if (not os.path.exists(tree_filepath)
			or not os.path.exists(sfreq_filepath)
			or not os.path.exists(dfreq_filepath)):
		print('Warning: Skipping \'' + filename_base + '\' as it is missing a tree or neighbourhood file.')
		return

	# get equilibrium frequencies
	with open(sfreq_filepath, 'r') as single_freq_file:
		single_frequencies = single_freq_file.readline()[:-1]
	with open(dfreq_filepath, 'r') as doublet_freq_file:
		doublet_frequencies = doublet_freq_file.readline()[:-1]

	for i in range(number):
		filename = filename_base + '_f' + str(i)
		length = random.randint(min_length, max_length)

		# generate ancestral sequence and secondary structure & convert to ct
		an_seq, ss_dbn = generate_sequence_structure_pair(length)
		ct = db_to_ct(ss_dbn, an_seq)

		# save generated sequence and structure (as dbn and ct)
		seq_out_filepath = os.path.join(outpath, 'sequences', filename + '.seq')
		with open(seq_out_filepath, 'w') as file:
			file.write(an_seq + '\n')
		with open(os.path.join(outpath, 'neighbourhoods/dbn', filename + '.dbn'), 'w') as file:
			file.write(an_seq + '\n' + ss_dbn + ' (0.0)\n')
		ct_out_filepath = os.path.join(outpath, 'neighbourhoods/ct', filename + '.ct')
		with open(ct_out_filepath, 'w') as file:
			file.write(ct + '\n')

		# generate alignment
		command = (sissi_filepath +
				   ' -k ' + str(seq_out_filepath) +
				   ' -nr' + str(ct_out_filepath) +
				   ' -fs ' + single_frequencies +
				   ' -fd ' + doublet_frequencies +
				   ' -l' + str(length) +
				   ' -a1' +
				   ' -oc' +
				   ' ' + str(tree_filepath))
		result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if len(result.stdout) == 0:
			print('Warning: Couldn\'t generate alignments for \'' + filename + '\':\n' + result.stderr)
			return

		# save alignment
		with open(os.path.join(outpath, 'alignments', filename + '.aln'), 'w') as outfile:
			outfile.write(result.stdout)

	print('Successfully generated famil' + ('ies' if number > 1 else 'y') + ' for \'' + filename_base + '\'.')


def generate_family_set(sissi_filepath, number, min_length, max_length, tree_dirpath, sfreq_dirpath, dfreq_dirpath, outpath):
	"""
	Generates <number> RNA families of uniformaly distributed lengths up to <maxlength> for each tree file in the given
	tree-directory, searching in the respectively given single- & doublet-frequencies-directories for files of the same
	name to use.
	For more information, refer to the generate_families() function.

		Parameters:
			sissi_filepath (str): Path to the compiled sissi099 file
			number (int): The number of families to generate
			min_length (int): Minimum allowed length of the ancestral sequences used to generate the families
			max_length (int): Maximum allowed length of the ancestral sequences used to generate the families
			tree_dirpath (str): Path to a directory containing tree files in the newick string format ('.seed_tree')
			sfreq_dirpath (str): Path to a directory containing files storing a single frequency vector ('.sfreq')
			dfreq_dirpath (str): Path to a directory containing files storing a doublet frequency vector ('.dfreq')
			outpath (str): Path to which to write the generated families
	"""

	# check paths
	if not os.path.exists(sissi_filepath):
		raise FileNotFoundError('No sissi099 file found at \'' + sissi_filepath + '\'')
	for path in [tree_dirpath, sfreq_dirpath, dfreq_dirpath]:
		if not os.path.exists(path):
			raise FileNotFoundError('No path named \'' + path + '\'')

	# make necessary directories
	os.makedirs(os.path.join(outpath, 'sequences'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/dbn'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/ct'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'alignments'), exist_ok=True)

	# generate
	filenames = os.listdir(tree_dirpath)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		generate_families(sissi_filepath,
						  number, min_length, max_length,
						  os.path.join(tree_dirpath, filename_base + '.seed_tree'),
						  os.path.join(sfreq_dirpath, filename_base + '.sfreq'),
						  os.path.join(dfreq_dirpath, filename_base + '.dfreq'),
						  outpath)
	print('Done.')


def setup_args(parser):
	parser.add_argument('--version', action='version', version='%(prog)s ' + '1.0')
	parser.add_argument('-sp', '--sissi-path', action='store', metavar='<str>',
						required=True, default=None,
						help='Provide the path to the compiled sissi099 file.')
	parser.add_argument('-t', '--type', metavar='<str>',
						required=True, default=None,
						choices=('f', 'families', 'a', 'alignments'),
						help='Select whether to generate whole families (\'f\', \'families\') or only alignments '
							 '(\'a\', \'alignments\')')
	parser.add_argument('-n', '--number', metavar='<int>',
						required=True, default=None,
						help='Provide the number of data points to generate.')
	parser.add_argument('-tp', '--tree-path', action='store', metavar='<str>',
						required=True, default=None,
						help='Provide the path to the directory containing newick string tree files (.seed_tree).')
	parser.add_argument('-sfp', '--single-freq-path', action='store', metavar='<str>',
						required=True, default=None,
						help='Provide the path to the directory containing single frequency vector files (.sfreq).')
	parser.add_argument('-dfp', '--doublet-freq-path', action='store', metavar='<str>',
						required=True, default=None,
						help='Provide the path to the directory containing doublet frequency vector files (.dfreq).')
	parser.add_argument('-op', '--out-path', action='store', metavar='<str>',
						required=True, default=None,
						help='Provide the path to the directory in which to save the generated data.')
	parser.add_argument('-fl', '--min-length', '--from-length', action='store', metavar='<int>',
						required=False, default=None,
						help='If generating families, provide the minimum allowed length for the generated sequences '
							 'and alignments')
	parser.add_argument('-tl', '--max-length', '--to-length', action='store', metavar='<int>',
						required=False, default=None,
						help='If generating families, provide the maximum allowed length for the generated sequences '
							 'and alignments')
	parser.add_argument('-np', '--neigh-path', action='store', metavar='<str>',
						required=False, default=None,
						help='If generating alignments, Provide the path to the directory containing dot bracket '
							 'notation secondary structure files (.dbn).')
	parser.add_argument('-ap', '--ali-path', action='store', metavar='<str>',
						required=False, default=None,
						help='If wanting to readd indels to generated alignments, provide the path to the directory '
							 'containing CLUSTAL format alignment files (.aln).')


def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	setup_args(parser)
	args = parser.parse_args()
	if args.type in ['a', 'alignments']:
		assert args.neigh_path is not None
		generate_alignment_set(args.sissi_path,
							   int(args.number),
							   args.tree_path,
							   args.neigh_path,
							   args.single_freq_path,
							   args.doublet_freq_path,
							   args.ali_path,
							   args.out_path)
	elif args.type in ['f', 'families']:
		assert args.min_length is not None
		assert args.max_length is not None
		generate_family_set(args.sissi_path,
							int(args.number),
							int(args.min_length),
							int(args.max_length),
							args.tree_path,
							args.single_freq_path,
							args.doublet_freq_path,
							args.out_path)
	else:
		raise argparse.ArgumentTypeError(
			'\'--type\' field required (options are \'a\', \'alignments\', \'f\', \'families\').')


if __name__ == '__main__':
	main()
