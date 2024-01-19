import os


def convert_fa_to_alidbn(infile, outpath):
	alioutpath = outpath + 'ali/'
	dbnoutpath = outpath + 'dbn/'

	try:
		os.makedirs(alioutpath)
	except FileExistsError:
		pass
	try:
		os.makedirs(dbnoutpath)
	except FileExistsError:
		pass

	ids = list()
	seqs = list()
	dbns = list()
	with open(infile) as seqdbnfile:
		line = seqdbnfile.readline()
		while line != '':
			if line[0] == '>':
				ids.append(line[:-1])
			elif line[0] in ['.', '(']:
				dbns.append(line[:-1])
			elif line[0] in ['A', 'C', 'G', 'U']:
				seqs.append(line[:-1])

			line = seqdbnfile.readline()

	for i in range(len(ids)):
		id = ids[i]
		if os.path.exists(alioutpath + id + '.ali'):
			id = id + '_2'
		seq = seqs[i]
		dbn = dbns[i]
		with open(alioutpath + id + '.ali', 'w') as outfile:
			outfile.write(' 1 ' + str(len(seq)) + '\n')
			outfile.write('_ ' + seq + '\n')
		with open(dbnoutpath + id + '.dbn', 'w') as outfile:
			outfile.write(seq + '\n')
			outfile.write(dbn + ' (0)\n')


def main():
	convert_fa_to_alidbn('uniform_len25-30_n10000.fa-train', 'seq_as_ali/')
	convert_fa_to_alidbn('uniform_len25-30_n2000.fa-test', 'seq_as_ali/')


if __name__ == '__main__':
	main()
