import os
import subprocess
import random


temp_name = ''


def generate_nucleotide(ambiguous_char):
	match ambiguous_char:
		case 'A' | 'C' | 'G' | 'U' | '-':
			return ambiguous_char
		case 'Y':
			return ['C', 'U'][random.randint(0, 1)]
		case 'R':
			return ['A', 'G'][random.randint(0, 1)]
		case 'W':
			return ['A', 'U'][random.randint(0, 1)]
		case 'S':
			return ['C', 'G'][random.randint(0, 1)]
		case 'K':
			return ['G', 'U'][random.randint(0, 1)]
		case 'M':
			return ['A', 'C'][random.randint(0, 1)]
		case 'D':
			return ['A', 'G', 'U'][random.randint(0, 2)]
		case 'V':
			return ['A', 'C', 'G'][random.randint(0, 2)]
		case 'H':
			return ['A', 'C', 'U'][random.randint(0, 2)]
		case 'B':
			return ['C', 'G', 'U'][random.randint(0, 2)]
		case 'N':
			return ['A', 'C', 'G', 'U'][random.randint(0, 3)]

	raise ValueError('\'' + ambiguous_char + '\' is not part of the IUPAC ambiguity code.')


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
				# TODO: remove
				if temp_name != '' and ali_name != temp_name:
					line = file.readline()
					continue
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


def stockholm_to_frequencies(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	# stockholm to ali_amb
	try:
		os.makedirs(outpath + 'ali_amb/')
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
				# TODO: remove
				if temp_name != '' and ali_name != temp_name:
					line = file.readline()
					continue
				with open(outpath + 'ali_amb/' + ali_name + '.ali', 'w') as outfile:
					number = len(ali)
					length = 0 if number == 0 else len(ali[0])
					outfile.write(' ' + str(number) + ' ' + str(length) + '\n')
					for seq in ali:
						outfile.write('_ ' + seq + '\n')

				ali = list()

			else:
				ali.append(line.split()[1])

			line = file.readline()

	# ali_amb to ali
	try:
		os.makedirs(outpath + 'ali/')
	except FileExistsError:
		pass
	filenames = os.listdir(outpath + 'ali_amb/')
	for filename in filenames:
		with open(outpath + 'ali_amb/' + filename, 'r') as file:
			with open(outpath + 'ali/' + filename, 'w') as outfile:
				line = file.readline()
				while line != '':
					split_line = line.split()
					if all([seg.isnumeric() for seg in split_line]):
						outfile.write(line)
					else:
						outfile.write(split_line[0] + ' ' + ''.join([generate_nucleotide(c) for c in split_line[1]]) + '\n')
					line = file.readline()

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
			while line != '':
				amounts[0] += line.count('A')
				amounts[1] += line.count('C')
				amounts[2] += line.count('G')
				amounts[3] += line.count('U')
				line = file.readline()
		count = sum(amounts)
		amounts = [amount / count for amount in amounts]
		with open(outpath + 'freq/' + filename.split('.')[0] + '.freq', 'w') as outfile:
			outfile.write(str(amounts[0] / count) + '  ' +
						  str(amounts[1] / count) + '  ' +
						  str(amounts[2] / count) + '  ' +
						  str(amounts[3] / count))


def fix_newick_string(filepath, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

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
				if char == ':' and fixed_newick_string_1[i+1].isdigit() and fixed_newick_string_1[i+2] == '.':  # TODO:?
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


def convert_rfam_data(seed_filepath, tree_path, tree_outpath, freq_outpath, nei_outpath):
	filenames = os.listdir(tree_path)
	for filename in filenames:
		# TODO: remove
		if temp_name != '' and filename.split('.')[0] != temp_name:
			continue
		fix_newick_string(tree_path + filename, tree_outpath)
	stockholm_to_frequencies(seed_filepath, freq_outpath)
	stockholm_to_neighbourhoods(seed_filepath, nei_outpath)
	print('Done.')


def main():
	convert_rfam_data('./data/rfam/Rfam_fixed.seed',
					  './data/rfam/seed_tree/original/',
					  './data/rfam/seed_tree/fixed/',
					  './data/rfam/seed_frequency/',
					  './data/rfam/seed_neighbourhood/')


if __name__ == '__main__':
	main()
