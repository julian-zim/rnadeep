import random
import os
import subprocess
import RNA


def db_to_ct(filepath, outpath):
	os.makedirs(outpath, exist_ok=True)

	with open(filepath, 'r') as file:
		seq = file.readline().split()[0]
		dbrs = file.readline().split()

	ptable = list(RNA.ptable(dbrs[0]))

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


def generate_ancestral_sequence(length=50):
	bases = ['A', 'C', 'G', 'U']
	seq = ''.join(random.choice(bases) for _ in range(length))
	return 'CAGUGCGAGCUA'
	#return seq
	pass


def generate_secondary_structure(sequence_filepath):
	with open(sequence_filepath, 'r') as file:
		seq = file.read().split()[0]
	dbrs, _ = RNA.fold(seq)
	return dbrs, seq


# sissi_filepath, tree_filepath, sfreq_dfilepath, dfreq_filepath, ali_filepath, outpath
def generate_family(out_filepath, length, as_filepath, ss_filepath, tree_filepath, sf_filepath, df_filepath, n):
	"""
	Generates n alignments using sissi for given equilibrium frequencies, neighbourhood system and phylogenetic tree.
	The raw alignments are used to re-add indels.

	Parameters:
	n (int): The number of alingments to generate
	directory (str): The path to a directory containing the information mentioned above in this exact folder structure:
		- seed_alignments: Directory with .aln Files, containing alignments
		- seed_frequencies:
			- single: Directory with .freq Files, containing unpaired nucleotide equilibrium frequencies for generation
			- doublet: Directory with .freq Files, containing paired nucleotide equilibrium frequencies for generation
		- seed_neighbourhoods:
			- nei: Directory with .nei files, containing the secondary consensus structure for generation
		- seed_trees:
			- rescaled: Directory with .seed_tree files, containing the phylogenetic tree for generation
	outpath (str): The path to which to write the generated alignments
	alignments (list): List of alignment file names to pick from the 'seed_alignments' folder for generation.
		None means all are picked.
	"""
	with open(as_filepath, 'w') as file:
		file.write(generate_ancestral_sequence(length) + '\n')

	ss = generate_secondary_structure(as_filepath)
	with open(ss_filepath, 'w') as file:
		file.write(ss[1] + '\n' + ss[0] + ' (0)\n')

	# secondary structure
	db_to_ct(ss_filepath, os.path.dirname(ss_filepath))
	ssct_filepath = os.path.join(os.path.dirname(ss_filepath), os.path.basename(ss_filepath).split('.')[0] + '.ct')

	# freqs
	with open(sf_filepath) as single_freq_file:
		single_frequencies = single_freq_file.readline()[:-1]
	with open(df_filepath) as doublet_freq_file:
		doublet_frequencies = doublet_freq_file.readline()[:-1]

	command = ('../../sissi099' +
			   ' -k ' + str(as_filepath) +
			   ' -nr' + str(ssct_filepath) +
			   ' -fs ' + single_frequencies +
			   ' -fd ' + doublet_frequencies +
			   ' -l' + str(length) +
			   ' -a' + str(n) +
			   ' -oc' +
			   ' ' + str(tree_filepath))
	result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if len(result.stdout) == 0:
		print('SISSI Error:\n' + result.stderr)
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

		for idx, ali in enumerate(ali_sets):
			with open(os.path.join(out_filepath + '_' + str(idx) + '.aln'), 'w') as outfile:
				outfile.write('CLUSTAL \n')
				for line in ali:
					outfile.write(line + '\n')
		print('Successfully generated alignment' + ('s' if n > 1 else '') + '.')


def main():
	filedir = 'Test/'
	filename = 'Test'
	generate_alignment('Test',
					   12,
					   filedir + filename + '.seq',
					   filedir + filename + '.dbn',
					   filedir + filename + '.seed_tree',
					   filedir + filename + '.sfreq',
					   filedir + filename + '.dfreq',
					   5)


if __name__ == '__main__':
	main()
