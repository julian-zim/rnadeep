import sys
import os
import subprocess
import RNA


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


def generate_alignments(sissi_filepath, n, tree_filepath, neigh_filepath, sfreq_dfilepath, dfreq_filepath, ali_filepath, outpath):
	"""
	Generates n RNA alignments using sissi for given equilibrium frequencies, neighbourhood system and phylogenetic tree.
	The raw alignments are used to re-add indels.
	Note: The given same neighbourhood will also be copied once for each generated alignment to create pairs for
	easier parsing into the network.

		Parameters:
			sissi_filepath (str): path to the compiled sissi099 file
			n (int): The number of alignments to generate
			tree_filepath (str): path to a tree file in the newick string format ('.seed_tree')
			neigh_filepath (str): path to a neighbourhood file in the sissi01 format ('.nei')
			sfreq_dfilepath (str): path to a file containing a single frequency vector ('.sfreq')
			dfreq_filepath (str): path to a file containing a doublet frequency vector ('.dfreq')
			ali_filepath (str): path to a file containing an alignment in the clustal format ('.aln')
			outpath (str): The path to which to write the generated alignments
	"""

	filename = os.path.basename(tree_filepath).split('.')[0]
	if (not os.path.exists(tree_filepath)
			or not os.path.exists(neigh_filepath)
			or not os.path.exists(sfreq_dfilepath)
			or not os.path.exists(dfreq_filepath)
			or not os.path.exists(ali_filepath)):
		print('Warning: Skipping \'' + filename + '\' as it is missing a tree, neighbourhood, frequency or alignment file.')
		return

	with open(sfreq_dfilepath, 'r') as single_freq_file:
		single_frequencies = single_freq_file.readline()[:-1]
	with open(dfreq_filepath, 'r') as doublet_freq_file:
		doublet_frequencies = doublet_freq_file.readline()[:-1]

	with open(neigh_filepath, 'r') as dbn_file:
		seq = dbn_file.readline()
		dbn = dbn_file.readline().split()[0]
	ct = db_to_ct(dbn, seq)

	for i in range(n):
		with open(os.path.join(outpath, 'neighbourhoods/dbn', filename + '_a' + str(i) + '.dbn'), 'w') as file:
			file.write(seq + dbn + ' (0.0)\n')
		with open(os.path.join(outpath, 'neighbourhoods/ct', filename + '_a' + str(i) + '.ct'), 'w') as file:
			file.write(ct + '\n')

	command = (sissi_filepath +
			   ' -nr' + str(neigh_filepath) +
			   ' -fs ' + single_frequencies +
			   ' -fd ' + doublet_frequencies +
			   ' -l' + str(len(dbn)) +
			   ' -a' + str(n) +
			   ' -oc' +
			   ' ' + str(tree_filepath))
	result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if len(result.stdout) == 0:
		print('Error in \'' + filename + '\':\n' + result.stderr)
	else:
		result_lines = result.stdout.split('\n')[:-1]  # for some reason a line without content is added in the end
		ali_sets = list()
		ali_idx = -1
		for line in result_lines:
			if line == 'CLUSTAL ':
				ali_sets.append(list())
				ali_idx += 1
			else:
				ali_sets[ali_idx].append(line)
		with open(ali_filepath, 'r') as alifile:
			alifile.readline()
			line = alifile.readline()
			line_idx = 0
			while line != '':
				for i, char in enumerate(line.split()[1]):
					if char == '-':
						for ali in ali_sets:
							seq = list(ali[line_idx])
							seq[i + len(ali[line_idx].split()[0]) + 1] = '-'
							ali[line_idx] = ''.join(seq)

				line = alifile.readline()
				line_idx += 1

		for idx, ali in enumerate(ali_sets):
			with open(os.path.join(outpath, 'alignments', filename + '_a' + str(idx) + '.aln'), 'w') as outfile:
				outfile.write('CLUSTAL \n')
				for line in ali:
					outfile.write(line + '\n')

		print('Successfully generated alignment' + ('s' if n > 1 else '') + ' for \'' + filename + '\'.')


def generate_alignment_set(sissi_filepath, n, tree_dirpath, neigh_dirpath, sfreq_dirpath, dfreq_dirpath, ali_dirpath, outpath):
	paths = [tree_dirpath, neigh_dirpath, sfreq_dirpath, dfreq_dirpath, ali_dirpath]
	for path in paths:
		if not os.path.exists(path):
			raise FileNotFoundError('No path named \'' + path + '\'')

	os.makedirs(os.path.join(outpath, 'alignments'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/ct'), exist_ok=True)
	os.makedirs(os.path.join(outpath, 'neighbourhoods/dbn'), exist_ok=True)

	filenames = os.listdir(tree_dirpath)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		generate_alignments(sissi_filepath,
							n,
							os.path.join(tree_dirpath, filename_base + '.seed_tree'),
							os.path.join(neigh_dirpath, filename_base + '.dbn'),
							os.path.join(sfreq_dirpath, filename_base + '.sfreq'),
							os.path.join(dfreq_dirpath, filename_base + '.dfreq'),
							os.path.join(ali_dirpath, filename_base + '.aln'),
							os.path.join(outpath))
	print('Done.')


def get_paths(rfam_path):
	tree_dirpath = os.path.join(rfam_path, 'seed_trees', 'rescaled')
	neigh_dirpath = os.path.join(rfam_path, 'seed_neighbourhoods', 'dbn')
	sfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'single')
	dfreq_dirpath = os.path.join(rfam_path, 'seed_frequencies', 'doublet')
	ali_dirpath = os.path.join(rfam_path, 'seed_alignments')
	return [tree_dirpath, neigh_dirpath, sfreq_dirpath, dfreq_dirpath, ali_dirpath]


def main():
	if len(sys.argv) != 9 and len(sys.argv) != 5:
		print('Usage: ./alignment_generator.py <sissi path> <amount> <trees path> <neighbourhoods path> <single frequencies path> <doublet frequencies path> <alignments path> <outpath>'
			  '\n       OR\n       ./alignment_generator.py <sissi path> <amount> <rfam path> <outpath>')
		return -1

	sissi_filepath = sys.argv[1]
	n = int(sys.argv[2])
	if len(sys.argv) == 5:
		paths = get_paths(sys.argv[3])
		outpath = sys.argv[4]
	elif len(sys.argv) == 9:
		paths = sys.argv[3:8]
		outpath = sys.argv[8]
	else:
		raise IOError('Impossible!')

	print('========== GENERATING ALIGNMENTS ==========')
	generate_alignment_set(sissi_filepath, n, *paths, outpath)


if __name__ == '__main__':
	main()
