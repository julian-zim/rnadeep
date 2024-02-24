import random
import sys
import os
import subprocess
import RNA


default_min_paired_sites = 25  # in percent


def db_to_ct(dbn, seq):
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


def generate_sequence_structure_pair(length=85, min_paired_sites=0):
	bases = ['A', 'C', 'G', 'U']
	seq = ''
	dbrs = ''
	paired_sites = -1
	while paired_sites < min_paired_sites:
		seq = ''.join(random.choice(bases) for _ in range(length))
		dbrs, _ = RNA.fold(seq)
		paired_sites = 0
		for char in dbrs:
			if char == '(' or char == ')':
				paired_sites += 1
	return seq, dbrs


def generate_family(sissi_filepath, n, length, tree_filepath, sfreq_filepath, dfreq_filepath, outpath):
	"""
	Generates n RNA families (consisting of an alignment and a secondary structure)
	for given equilibrium frequencies and a phylogenetic tree, using:
	- a random ancestral sequence
	- RNAfold to predict a consensus structure
	- SISSI simulate a corresponding homologous sequence alignment.

	Parameters:
	sissi_filepath (str): path to the compiled sissi099 file
	n (int): The number of families to generate
	length(int): Length of the ancestral sequence used to generate the family
	tree_filepath (str): path to a tree file in the newick string format ('.seed_tree')
	sfreq_dfilepath (str): path to a file containing a single frequency vector ('.sfreq')
	dfreq_filepath (str): path to a file containing a doublet frequency vector ('.dfreq')
	outpath (str): The path to which to write the generated families
	"""

	filename_base = os.path.basename(tree_filepath).split('.')[0]
	if (not os.path.exists(tree_filepath)
			or not os.path.exists(sfreq_filepath)
			or not os.path.exists(dfreq_filepath)):
		print('Warning: Skipping \'' + filename_base + '\' as it is missing a tree or neighbourhood file.')
		return

	with open(sfreq_filepath, 'r') as single_freq_file:
		single_frequencies = single_freq_file.readline()[:-1]
	with open(dfreq_filepath, 'r') as doublet_freq_file:
		doublet_frequencies = doublet_freq_file.readline()[:-1]

	for i in range(n):
		filename = filename_base + '_f' + str(i)

		an_seq, ss_dbn = generate_sequence_structure_pair(length, int(length * default_min_paired_sites / 100))
		ct = db_to_ct(ss_dbn, an_seq)

		seq_out_filepath = os.path.join(outpath, 'sequences', filename + '.seq')
		with open(seq_out_filepath, 'w') as file:
			file.write(an_seq + '\n')
		with open(os.path.join(outpath, 'neighbourhoods/dbn', filename + '.dbn'), 'w') as file:
			file.write(an_seq + '\n' + ss_dbn + ' (0.0)\n')
		ct_out_filepath = os.path.join(outpath, 'neighbourhoods/ct', filename + '.ct')
		with open(ct_out_filepath, 'w') as file:
			file.write(ct + '\n')

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
			print('SISSI Error:\n' + result.stderr)
		else:
			with open(os.path.join(outpath, 'alignments', filename + '.aln'), 'w') as outfile:
				outfile.write(result.stdout)

	print('Successfully generated famil' + ('ies' if n > 1 else 'y') + '.')


def generate_family_set(sissi_filepath, n, length, tree_dirpath, sfreq_dirpath, dfreq_dirpath, outpath):
	paths = [tree_dirpath, sfreq_dirpath, dfreq_dirpath]
	for path in paths:
		if not os.path.exists(path):
			raise FileNotFoundError('No path named \'' + path + '\'')

	os.makedirs(os.path.join(outpath, 'sequences'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'alignments'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/dbn'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/ct'), exist_ok=True)

	filenames = os.listdir(tree_dirpath)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		generate_family(sissi_filepath,
						n,
						length,
						os.path.join(tree_dirpath, filename_base + '.seed_tree'),
						os.path.join(sfreq_dirpath, filename_base + '.sfreq'),
						os.path.join(dfreq_dirpath, filename_base + '.dfreq'),
						outpath)
	print('Done.')


def get_paths(rfam_path):
	tree_dirpath = os.path.join(rfam_path, 'seed_trees', 'rescaled')
	sfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'single')
	dfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'doublet')
	return [tree_dirpath, sfreq_dirpath, dfreq_dirpath]


def main():
	if len(sys.argv) != 8 and len(sys.argv) != 6:
		print('Usage: ./alignment_generator.py <sissi path> <amount> <length> <trees path> <single frequencies path> <doublet frequencies path> <outpath>'
			  '\n       OR\n       ./alignment_generator.py <sissi path> <amount> <length> <rfam path> <outpath>')
		return -1

	sissi_filepath = sys.argv[1]
	n = int(sys.argv[2])
	length = int(sys.argv[3])
	if len(sys.argv) == 6:
		paths = get_paths(sys.argv[4])
		outpath = sys.argv[5]
	elif len(sys.argv) == 8:
		paths = sys.argv[4:7]
		outpath = sys.argv[7]
	else:
		raise IOError('Impossible!')

	print('========== GENERATING FAMILIES ==========')
	generate_family_set(sissi_filepath, n, length, *paths, outpath)


if __name__ == '__main__':
	main()
