import os
import subprocess


def ct_to_nei(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	with open(inpath + filename, 'r') as file:
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


def dotbracket_to_ct(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	result = subprocess.run('RNAfold --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		raise RuntimeError('ViennaRNA is not installed.')

	command = 'b2ct < ' + inpath + filename + '> ' + outpath + filename + '.ct'
	result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		print('Warning: Couldn\'t convert file ' + filename)
		os.remove(outpath + filename + '.ct')


def wuss_to_dotbracket(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	with open(inpath + filename, 'r') as file:
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

	with open(outpath + filename, 'w') as outfile:
		outfile.write(fixed_content)


def stockholm_to_wuss(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	with open(inpath + filename, 'r') as file:
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
				with open(outpath + ali_name + '.fa', 'w') as outfile:
					outfile.write(ali_cons_sequence + '\n')
					outfile.write(ali_cons_structure + '\n')

			else:
				pass  # actual alignment data. only important in stockholm_to_frequencies function, not here

			line = file.readline()


def stockholm_to_neighbourhoods(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	stockholm_to_wuss(inpath, filename, outpath + 'wuss/')

	filenames = os.listdir(outpath + 'wuss/')
	for filename in filenames:
		wuss_to_dotbracket(outpath + 'wuss/', filename, outpath + 'dotbracket/')

	filenames = os.listdir(outpath + 'dotbracket/')
	for filename in filenames:
		dotbracket_to_ct(outpath + 'dotbracket/', filename, outpath + 'ct/')

	filenames = os.listdir(outpath + 'ct/')
	for filename in filenames:
		ct_to_nei(outpath + 'ct/', filename, outpath + 'nei/')


def stockholm_to_frequencies(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	# stockholm to ali
	try:
		os.makedirs(outpath + 'fa/')
	except FileExistsError:
		pass
	with open(inpath + filename, 'r') as file:
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
				with open(outpath + 'fa/' + ali_name + '.fa', 'w') as outfile:
					for seq in ali:
						outfile.write(seq + '\n')

				ali = list()

			else:
				ali.append(line.split()[1])

			line = file.readline()

	# ali to freq
	try:
		os.makedirs(outpath + 'freq/')
	except FileExistsError:
		pass
	filenames = os.listdir(outpath + 'fa/')
	for filename in filenames:
		amounts = [0, 0, 0, 0]
		with open(outpath + 'fa/' + filename, 'r') as file:
			line = file.readline()
			while line != '':
				amounts[0] += line.count('A')
				amounts[1] += line.count('C')
				amounts[2] += line.count('G')
				amounts[3] += line.count('U')
				line = file.readline()
		count = sum(amounts)
		amounts = [amount / count for amount in amounts]
		with open(outpath + 'freq/' + filename + '.freq', 'w') as outfile:
			outfile.write(str(amounts[0] / count) + '  ' +
						  str(amounts[1] / count) + '  ' +
						  str(amounts[2] / count) + '  ' +
						  str(amounts[3] / count))


def fix_newick_string(inpath, filename, outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	with open(inpath + filename, 'r') as file:
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
		for char in fixed_newick_string_1:

			if predict_name:
				if char != '(' and char != ',':
					name = True
					predict_name = False

			if name:
				if char == '(' or char == ')':
					fixed_newick_string_2 += '_'  # replace brackets in names with underscores
				elif char == ':':
					name = False
					fixed_newick_string_2 += char
				else:
					fixed_newick_string_2 += char

			else:
				if not predict_name and (char == '(' or char == ','):
					predict_name = True
				fixed_newick_string_2 += char

	with open(outpath + filename, 'w') as outfile:
		outfile.write(fixed_newick_string_2)


def main():
	stockholm_to_frequencies('./data/', 'Rfam_fixed.seed', './data/Rfam.seed_frequency/')
	stockholm_to_neighbourhoods('./data/', 'Rfam_fixed.seed', './data/Rfam.seed_neighbourhood/')
	filenames = os.listdir('./data/Rfam.seed_tree/original/')
	for filename in filenames:
		# TODO: remove
		if temp_name != '' and filename.split('.')[0] != temp_name:
			continue
		fix_newick_string('./data/Rfam.seed_tree/original/', filename, './data/Rfam.seed_tree/fixed/')


if __name__ == '__main__':
	temp_name = ''
	main()
