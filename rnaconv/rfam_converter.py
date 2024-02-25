import sys
import os
import subprocess
import numpy as np
import textdistance
import shutil
import RNA
from ete3 import TreeNode


def rescale_newick_strings(treedirpath, alidirpath, outpath):
	"""
	Rescales the tree branch lengths for trees which corresponding sequence alignments sequences are over 95% similar
	with respect to their mean pairwise hamming distance, in order to increase the evolution rate when using the tree
	for evolutionary simulation.
	The rescale factor is 2.

		Parameters:
			treedirpath (str): path to the directory containing the tree files in newick string format
			alidirpath (str): path to the directory containing the alignment files in CLUSTAL format
			outpath (str): path to the directory in which to save the rescaled trees in the newick string format.
	"""

	os.makedirs(outpath, exist_ok=True)

	tree_filenames = os.listdir(treedirpath)
	for tree_filename in tree_filenames:
		filename = tree_filename.split('.')[0]

		with open(os.path.join(alidirpath, filename + '.aln'), 'r') as file:
			ali = [line.split()[1] for line in file.readlines()[1:]]

		pdists = np.zeros((len(ali), len(ali)))
		for i, seq1 in enumerate(ali):
			for j, seq2 in enumerate(ali):
				pdists[i, j] = textdistance.hamming(seq1, seq2)
		mean_pw_distance = np.mean(pdists)

		shutil.copy(os.path.join(treedirpath, tree_filename), outpath)
		if mean_pw_distance / len(ali[0]) < 0.05:
			shutil.copy(os.path.join(treedirpath, tree_filename), outpath)
			with open(os.path.join(outpath, tree_filename), 'r') as file:
				newick_string = file.read()
			tree = TreeNode(newick_string)
			for node in tree.traverse():
				node.dist *= 2  # rescale factor
			tree_rescaled = tree.write(format=1)
			with open(os.path.join(outpath, filename.split('.')[0] + '.seed_tree'), 'w+') as file:
				file.write(tree_rescaled)
			print('Info: Doubled branch lengths for tree \'' + filename + '\' due to the corresponding alignment\'s '
				  'sequences being 95% similar.')


def fix_newick_strings(treedirpath, outpath):
	"""
	Fixes newick strings by replacing every control character (e.g. '(', ')', ',', '.', ':') within a node name with an
	underscore.
	Additionally, multifurcations are resolved and non-leaf node labels are removed.

	These three steps are nessecary for SISSI to be able to parse the Rfam tree files.

		Parameters:
			treedirpath (str): path to the directory containing the tree files in newick string format
			outpath (str): path to the directory in which to save the trees in the fixed newick string format.
	"""

	os.makedirs(outpath, exist_ok=True)

	filenames = os.listdir(treedirpath)
	for filename in filenames:
		filepath = os.path.join(treedirpath, filename)

		with (open(filepath, 'r') as file):
			newick_string = file.read()

		# fix control characters in names
		fixed_newick_string_0 = ''
		predict_name = False
		name = False

		for i, char in enumerate(newick_string):

			if predict_name:
				if char != '(' and char != ',':
					name = True
					predict_name = False

			if name:
				if char == ':' and newick_string[i + 1].isdigit() and newick_string[i + 2] == '.':
					name = False
					fixed_newick_string_0 += char
				else:
					if char == '(' or char == ')' or char == ',' or char == ':' or char == '.':
						char_to_add = '_'  # replace brackets, [semi-]colons and kommas in names with underscores
					else:
						char_to_add = char
					fixed_newick_string_0 += char_to_add

			else:
				if not predict_name and (char == '(' or char == ','):
					predict_name = True
				fixed_newick_string_0 += char

		# resolve multifurcations & non-leaf node names
		tree = TreeNode(fixed_newick_string_0)
		tree.resolve_polytomy()
		fixed_newick_string_1 = tree.write(format=1)
		'''
		# remove non-leaf node names (for whatever reason SISSi cant deal with them)
		fixed_newick_string_2 = ''
		alert = False
		for char in fixed_newick_string_1:
			if alert:
				if char == ':' or char == ';':
					fixed_newick_string_2 += char
					alert = False
				else:
					pass  # remove non-leaf names
			else:
				if char == ')':
					alert = True
				fixed_newick_string_2 += char
		'''

		# save
		with open(os.path.join(outpath, filepath.split('/')[-1]), 'w') as outfile:
			outfile.write(fixed_newick_string_1)


def obtain_equilibrium_frequencies(alidirpath, neighdirpath, outpath):
	"""
	Extracts the equilibrium frequencies for unpaired single nucleotides and nucleotide pairs from an alignment.

	It counts the occurences of single nucleotides per unpaired site and saves them in a 4-vector.
	Then It counts the occurences of nucleotide pairs per paired site tuple and saves them in a 16-vector.
	Then it adds pseudocounts (+1 for each element) and normalizes.

		Parameters:
			alidirpath (str): path to the directory containing the alignment files in CLUSTAL format
			neighdirpath (str): path to the directory containing the alignment consensus structure files in sissi0.1
				(.nei) format
			outpath (str): path to the directory in which to save the extracted unpaired single and paired nucleotide
				equilibrium frequencies
	"""

	os.makedirs(outpath, exist_ok=True)
	os.makedirs(os.path.join(outpath, 'single'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'doublet'), exist_ok=True)

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
		'V': np.array([1./3., 1./3., 1./3., 0]),
		'H': np.array([1./3., 1./3., 0, 1./3.]),
		'D': np.array([1./3., 0, 1./3., 1./3.]),
		'B': np.array([0, 1./3., 1./3., 1./3.]),
		'N': np.array([0.25, 0.25, 0.25, 0.25]),
	}

	for alifilename in os.listdir(alidirpath):
		filename_base = alifilename.split('.')[0]
		if not os.path.exists(os.path.join(neighdirpath, 'nei', filename_base + '.nei')):
			print('Warning: Couldn\'t obtain equilibrium frequencies of \"' + filename_base +
				  '\" since there is no neighbourhood information.')
			continue

		alignment = list()
		with open(os.path.join(alidirpath, alifilename)) as alifile:
			alifile.readline()
			line = alifile.readline()
			while line != '':
				alignment.append(line.split()[1])
				line = alifile.readline()

		neighbourhood = set()
		with open(os.path.join(neighdirpath, 'nei', filename_base + '.nei')) as neifile:
			line = neifile.readline()
			while line != '':
				split_line = line.split()
				if len(split_line) < 3:
					split_line.append('-1')
				pair = tuple(sorted((int(split_line[1].split('|')[0]), int(split_line[2].split(':')[0]))))
				neighbourhood.add(pair)
				line = neifile.readline()

		single_frequencies = np.ones(4)  # ones because of pseudocounts
		doublet_frequencies = np.ones(16)  # ones because of pseudocounts
		for pair in neighbourhood:
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

		# normalization
		sum_single_freq = sum(single_frequencies)
		if sum_single_freq > 0:
			single_frequencies /= sum_single_freq
		sum_double_freq = sum(doublet_frequencies)
		if sum_double_freq > 0:
			doublet_frequencies /= sum_double_freq

		with open(os.path.join(outpath, 'single', filename_base + '.sfreq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in single_frequencies]) + '\n')
		with open(os.path.join(outpath, 'doublet', filename_base + '.dfreq'), 'w') as outfile:
			outfile.write(' '.join([str(e) for e in doublet_frequencies]) + '\n')


def ct_to_nei(filepath, outpath):
	"""
	Converts the consensus structures contained in the connect table input file into the sissi0.1 (.nei) format
	Note: SISSI can also use connect table files directly, but this filetype is used for obtaining the equilibrium
	frequencies later on.

		Parameters:
			filepath (str): path to the file in connect table format, containing the secondary structure
			outpath (str): path to the directory in which to save the resulting sissi0.1 (.nei) file
	"""

	os.makedirs(outpath, exist_ok=True)

	filename = filepath.split('/')[-1].split('.')[0]

	with open(filepath, 'r') as file:
		with open(os.path.join(outpath, filename + '.nei'), 'w') as outfile:
			file.readline()
			line = file.readline()
			while line != '':
				data = line.split()
				pos1 = str(int(data[0]) - 1)
				pos2 = str(int(data[4]) - 1)
				gap = ' ' * (5 - len(str(pos1)))
				outfile.write('Pos' + gap + pos1 + '|' + ((' ' + pos2 + ':(1.000000)\n') if pos2 != '-1' else '\n'))
				line = file.readline()


def db_to_ct(filepath, outpath):
	"""
	Converts the consensus structures contained in the dot bracket notation input file into the connect table format

		Parameters:
			filepath (str): path to the file in dot bracket notation format, containing the secondary structure
			outpath (str): path to the directory in which to save the resulting connect table file
	"""

	os.makedirs(outpath, exist_ok=True)

	with open(filepath, 'r') as file:
		seq = file.readline().split()[0]
		dbrs = file.readline().split()

	ptable = list(RNA.ptable(dbrs[0]))
	'''result = subprocess.run('RNAfold --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		raise RuntimeError('ViennaRNA is not installed.')

	command = 'b2ct < ' + filepath + ' > ' + str(os.path.join(outpath, filename + '.ct'))  # TODO: use library
	result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		print('Warning: Couldn\'t convert file ' + filename + ' to ct. Traceback: ' + result.stderr.strip('\n'))
		os.remove(os.path.join(outpath, filename + '.ct'))'''

	column1 = list(range(1, ptable[0] + 1))
	column2 = [c for c in seq]
	column3 = [idx - 1 for idx in column1[:-1]] + [0]
	column4 = [idx + 1 for idx in column1]
	column5 = ptable[1:]
	column6 = column1
	with open(os.path.join(outpath, os.path.basename(filepath).split('/')[-1].split('.')[0] + '.ct'), 'w') as file:
		file.write('{: >5}'.format(ptable[0]) + ' ENERGY =     0.0    1\n')
		for i in column3:
			file.write('{: >5}'.format(column1[i]) +
					   ' ' + column2[i] +
					   '{: >8}'.format(column3[i]) +
					   '{: >5}'.format(column4[i]) +
					   '{: >5}'.format(column5[i]) +
					   '{: >5}'.format(column6[i]) + '\n')


def wuss_to_db(filepath, outpath):
	"""
	Converts the consensus structures contained in the wuss input file into the dot bracket notation format

		Parameters:
			filepath (str): path to the file in wuss format, containing the secondary structure
			outpath (str): path to the directory in which to save the resulting dot bracket notation file
	"""

	os.makedirs(outpath, exist_ok=True)

	filename = filepath.split('/')[-1].split('.')[0]

	with open(filepath, 'r') as file:
		content = file.read()
		content_lines = content.split('\n')

		fixed_line = RNA.db_from_WUSS(content_lines[1])
		'''fixed_line = ''
		for char in content_lines[1]:
			if char == '<' or char == '{' or char == '[':
				fixed_line += '('
			elif char == '>' or char == '}' or char == ']':
				fixed_line += ')'
			elif char == ',' or char == '-' or char == '_' or char == '~' or char == ':' \
				or char == 'A' or char == 'a' or char == 'B' or char == 'b' \
				or char == 'C' or char == 'c' or char == 'D' or char == 'd':
				fixed_line += '.'
			elif char == '(' or char == ')' or char == '.':
				fixed_line += char
			else:
				raise ValueError('WRONG FORMAT: ' + char + ', in ' + filename)'''

		fixed_content = content_lines[0] + '\n' + fixed_line + ' (0)\n'

	with open(os.path.join(outpath, filename + '.dbn'), 'w') as outfile:
		outfile.write(fixed_content)


def stockholm_to_wuss(filepath, outpath):
	"""
	Converts the consensus structures contained in the STOCKHOLM input file into single files in the wuss format

		Parameters:
			filepath (str): path to the Rfam.seed file in STOCKHOLM format, containing the families
			outpath (str): path to the directory in which to save the resulting wuss file
	"""

	os.makedirs(outpath, exist_ok=True)

	with open(filepath, 'r') as file:
		ali_name = ''
		ali_cons_structure = ''

		line = file.readline()
		while line != '':

			if line[0] == '\n':
				pass

			elif line[0] == '#':
				if line[1:4] == '=GF':
					split_line = line.split()
					if split_line[1] == 'AC':
						ali_name = split_line[2]
					else:
						pass  # skip other compulsory fields
				elif line[1:4] == '=GC':
					split_line = line.split()
					if split_line[1] == 'SS_cons':
						ali_cons_structure = split_line[2]
					elif split_line[1] == 'RF':
						ali_cons_sequence = split_line[2].upper()
					else:
						pass  # skip other compulsory fields

			elif line[0:2] == '//':
				with open(os.path.join(outpath, ali_name + '.wuss'), 'w') as outfile:
					outfile.write(ali_cons_sequence + '\n')
					outfile.write(ali_cons_structure + '\n')

			else:
				pass  # actual alignment data. only important in stockholm_to_frequencies function, not here

			line = file.readline()


def stockholm_to_neighbourhoods(filepath, outpath):
	"""
	Converts the consensus structures contained in the STOCKHOLM input file into single files in the following formats:
		- wuss
		- dbn
		- ct
		- nei

		Parameters:
			filepath (str): path to the Rfam.seed file in STOCKHOLM format, containing the families
			outpath (str): path to the directory in which to save the extracted consensus structures
	"""

	os.makedirs(outpath, exist_ok=True)

	stockholm_to_wuss(filepath, os.path.join(outpath, 'wuss'))

	filenames = os.listdir(os.path.join(outpath, 'wuss'))
	for filename in filenames:
		wuss_to_db(os.path.join(outpath, 'wuss', filename), os.path.join(outpath, 'dbn'))

	filenames = os.listdir(os.path.join(outpath, 'dbn'))
	for filename in filenames:
		db_to_ct(os.path.join(outpath, 'dbn', filename), os.path.join(outpath, 'ct'))

	filenames = os.listdir(os.path.join(outpath, 'ct'))
	for filename in filenames:
		ct_to_nei(os.path.join(outpath, 'ct', filename), os.path.join(outpath, 'nei'))


def stockholm_to_alignments(filepath, outpath):
	"""
	Converts the alignments contained in the STOCKHOLM input file into CLUSTAL files

		Parameters:
			filepath (str): path to the Rfam.seed file in STOCKHOLM format, containing the families
			outpath (str): path to the directory in which to save the extracted alignments in the CLUSTAL format
	"""

	os.makedirs(outpath, exist_ok=True)

	with open(filepath, 'r') as file:
		ali_name = ''
		ali = list()

		line = file.readline()
		while line != '':

			if line[0] == '\n':
				pass

			elif line[0] == '#':
				if line[1:4] == '=GF':
					split_line = line.split()
					if split_line[1] == 'AC':
						ali_name = split_line[2]
					else:
						pass  # skip other compulsory fields
				else:
					pass  # skip other markups

			elif line[0:2] == '//':
				with open(os.path.join(outpath, ali_name + '.aln'), 'w') as outfile:
					#number = len(ali)
					#length = 0 if number == 0 else len(ali[0])
					#outfile.write(' ' + str(number) + ' ' + str(length) + '\n')
					outfile.write('CLUSTAL ' + '\n')
					for seq in ali:
						outfile.write(seq + '\n')

				ali = list()

			else:
				split_line = line.split()
				ali.append(split_line[0] + ' ' + split_line[1])

			line = file.readline()


def convert_rfam_data(seed_filepath,
					  ali_outpath,
					  neigh_outpath,
					  freq_outpath,
					  tree_path,
					  tree_fixed_outpath,
					  tree_rescaled_outpath):
	"""
	Calls the nessecary functions to convert the whole rfam database into single files, preparing them to be used by
	SISSI.

		Parameters:
			seed_filepath (str): path to the Rfam.seed file in STOCKHOLM format, containing the families
			ali_outpath (str): path to the directory in which to save the extracted alignments in the CLUSTAL format
			neigh_outpath (str): path to the directory in which to save the extracted consensus structures in the
				- wuss
				- dbn
				- ct
				- nei
				formats, respectively
			freq_outpath (str): path to the directory in which to save the extracted paired and unpaired nucleotide
				equilibrium frequencies
			tree_path (str): path to the directory in which Rfam tree files in newick format (.seed_tree) files are
				located
			tree_fixed_outpath (str): path to the directory in which to save the fixed tree newick strings
			tree_rescaled_outpath (str): path to the directory in which to save the rescaled tree newick strings
	"""

	stockholm_to_alignments(seed_filepath, ali_outpath)
	stockholm_to_neighbourhoods(seed_filepath, neigh_outpath)
	obtain_equilibrium_frequencies(ali_outpath, neigh_outpath, freq_outpath)
	fix_newick_strings(tree_path, tree_fixed_outpath)
	rescale_newick_strings(tree_fixed_outpath, ali_outpath, tree_rescaled_outpath)
	print('Done.')


def main():
	if len(sys.argv) != 2:
		print('Usage: ./rfam_converter.py <rfam_path>')
		return -1

	rfam_path = sys.argv[1]

	if (not os.path.exists(os.path.join(rfam_path, 'Rfam.seed'))
			or not os.path.exists(os.path.join(rfam_path, 'seed_trees/original'))):
		print('Required rfam folder structure:\n- rfam\n\tRfam.seed\n\t- seed_trees\n\t\t- original\n\t\t\t<tree files>')
		return -1

	print('========== CONVERTING RFAM DATA POINTS ==========')
	convert_rfam_data(os.path.join(rfam_path, 'Rfam.seed'),
					  os.path.join(rfam_path, 'seed_alignments'),
					  os.path.join(rfam_path, 'seed_neighbourhoods'),
					  os.path.join(rfam_path, 'seed_frequencies'),
					  os.path.join(rfam_path, 'seed_trees/original'),
					  os.path.join(rfam_path, 'seed_trees/fixed'),
					  os.path.join(rfam_path, 'seed_trees/rescaled'))


if __name__ == '__main__':
	main()
