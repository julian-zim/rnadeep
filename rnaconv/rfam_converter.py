import os
import subprocess
import numpy as np
import math


"""
def obtain_single_frequencies(outpath):
	# ali to freq
	try:
		os.makedirs(outpath + 'freq/')
	except FileExistsError:
		pass
	filenames = os.listdir(outpath + 'ali/')
	for filename in filenames:
		amounts = [0, 0, 0, 0]
		with open(outpath + 'ali/' + filename, 'r') as file:
			line = file.readline()
			while line != '':  # TODO: add ambiguity amounts
				amounts[0] += line.count('A')
				amounts[1] += line.count('C')
				amounts[2] += line.count('G')
				amounts[3] += line.count('U')
				line = file.readline()
		count = sum(amounts)
		amounts = [amount / count for amount in amounts]
		with open(outpath + 'freq/' + filename.split('.')[0] + '.freq', 'w') as outfile:
			outfile.write(str(amounts[0]) + '  ' +
						  str(amounts[1]) + '  ' +
						  str(amounts[2]) + '  ' +
						  str(amounts[3]))
"""


def obtain_doublet_frequencies(alidirpath, neighdirpath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

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
		if not os.path.exists(neighdirpath + 'nei/' + filename_base + '.nei'):
			print('Warning: Couldn\'t obtain doublet frequencies of \"' + filename_base +
				  '\" since there is no neighbourhood information.')
			continue

		alignment = list()
		with open(alidirpath + alifilename) as alifile:
			alifile.readline()
			line = alifile.readline()
			while line != '':
				alignment.append(line.split()[1])
				line = alifile.readline()

		neighbourhood = set()
		with open(neighdirpath + 'nei/' + filename_base + '.nei') as neifile:
			line = neifile.readline()
			while line != '':
				split_line = line.split()
				if len(split_line) < 3:
					line = neifile.readline()
					continue
				pair = tuple(sorted((int(split_line[1].split('|')[0]), int(split_line[2].split(':')[0]))))
				neighbourhood.add(pair)
				line = neifile.readline()

		doublet_frequencies = np.zeros(16)
		for pair in neighbourhood:
			for seq in alignment:
				char_i, char_j = seq[pair[0]].upper(), seq[pair[1]].upper()
				if char_i != '-' and char_j != '-':
					summand = np.outer(base_to_ids[char_i], base_to_ids[char_j]).flatten()
					doublet_frequencies += summand
		sum_freq = sum(doublet_frequencies)
		if sum_freq > 0:
			doublet_frequencies /= sum_freq

		with open(outpath + filename_base + '.freq', 'w') as outfile:
			outfile.write(' '.join([str(e) for e in doublet_frequencies]) + '\n')


def ct_to_nei(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	filename = filepath.split('/')[-1].split('.')[0]

	with open(filepath, 'r') as file:
		with open(outpath + filename + '.nei', 'w') as outfile:
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
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	filename = filepath.split('/')[-1].split('.')[0]

	result = subprocess.run('RNAfold --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		raise RuntimeError('ViennaRNA is not installed.')

	command = 'b2ct < ' + filepath + '> ' + outpath + filename + '.ct'
	result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		print('Warning: Couldn\'t convert file ' + filename + ' to ct.')
		os.remove(outpath + filename + '.ct')


def wuss_to_db(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	filename = filepath.split('/')[-1].split('.')[0]

	with open(filepath, 'r') as file:
		content = file.read()

		content_lines = content.split('\n')
		fixed_line = ''
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
				raise ValueError('WRONG FORMAT: ' + char + ', in ' + filename)

		fixed_content = content_lines[0] + '\n' + fixed_line + ' (0)\n'

	with open(outpath + filename + '.dbn', 'w') as outfile:
		outfile.write(fixed_content)


def stockholm_to_wuss(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

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
				with open(outpath + ali_name + '.wuss', 'w') as outfile:
					outfile.write(ali_cons_sequence + '\n')
					outfile.write(ali_cons_structure + '\n')

			else:
				pass  # actual alignment data. only important in stockholm_to_frequencies function, not here

			line = file.readline()


def stockholm_to_neighbourhoods(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	stockholm_to_wuss(filepath, outpath + 'wuss/')

	filenames = os.listdir(outpath + 'wuss/')
	for filename in filenames:
		wuss_to_db(outpath + 'wuss/' + filename, outpath + 'dbn/')

	filenames = os.listdir(outpath + 'dbn/')
	for filename in filenames:
		db_to_ct(outpath + 'dbn/' + filename, outpath + 'ct/')

	filenames = os.listdir(outpath + 'ct/')
	for filename in filenames:
		ct_to_nei(outpath + 'ct/' + filename, outpath + 'nei/')


def stockholm_to_alignments(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass
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
				with open(outpath + ali_name + '.ali', 'w') as outfile:
					number = len(ali)
					length = 0 if number == 0 else len(ali[0])
					outfile.write(' ' + str(number) + ' ' + str(length) + '\n')
					for seq in ali:
						outfile.write('_ ' + seq + '\n')

				ali = list()

			else:
				ali.append(line.split()[1])

			line = file.readline()


def fix_newick_strings(dirpath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	filenames = os.listdir(dirpath)
	for filename in filenames:
		filepath = dirpath + filename

		with (open(filepath, 'r') as file):
			newick_string = file.read()

			fixed_newick_string_1 = ''
			alert = False
			for char in newick_string:
				if alert:
					if char == ':' or char == ';':
						fixed_newick_string_1 += char
						alert = False
					else:
						pass  # remove non-leaf names
				else:
					if char == ')':
						alert = True
					fixed_newick_string_1 += char

			fixed_newick_string_2 = ''
			predict_name = False
			name = False
			for i, char in enumerate(fixed_newick_string_1):

				if predict_name:
					if char != '(' and char != ',':
						name = True
						predict_name = False

				if name:
					if char == ':' and fixed_newick_string_1[i + 1].isdigit() and fixed_newick_string_1[i + 2] == '.':
						name = False
						fixed_newick_string_2 += char
					elif char == '(' or char == ')' or char == ',' or char == ':':
						fixed_newick_string_2 += '_'  # replace brackets and colons in names with underscores
					else:
						fixed_newick_string_2 += char

				else:
					if not predict_name and (char == '(' or char == ','):
						predict_name = True
					fixed_newick_string_2 += char

		with open(outpath + filepath.split('/')[-1], 'w') as outfile:
			outfile.write(fixed_newick_string_2)


def convert_rfam_data(tree_path, tree_outpath, seed_filepath, ali_outpath, neigh_outpath, freq_outpath):
	fix_newick_strings(tree_path, tree_outpath)
	stockholm_to_alignments(seed_filepath, ali_outpath)
	stockholm_to_neighbourhoods(seed_filepath, neigh_outpath)
	obtain_doublet_frequencies(ali_outpath, neigh_outpath, freq_outpath)
	print('Done.')


def main():
	convert_rfam_data('./data/ambiguous/rfam/seed_trees/original/',
					  './data/ambiguous/rfam/seed_trees/fixed/',
					  './data/ambiguous/rfam/Rfam.seed',
					  './data/ambiguous/rfam/seed_alignments/',
					  './data/ambiguous/rfam/seed_neighbourhoods/',
					  './data/ambiguous/rfam/seed_frequencies/')


if __name__ == '__main__':
	main()
