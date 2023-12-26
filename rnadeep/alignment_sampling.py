import os


def parse_seed_alignment(ali_path, dbs_path, filename):
	with open(ali_path + filename + '.ali') as ali_file:
		alignment = list()
		line = ali_file.readline()
		while line != '':
			alignment.append(line[:-1])
			line = ali_file.readline()
	with open(dbs_path + filename + '.dbs') as dbs_file:
		dbs = dbs_file.readlines()[1].split()[0]

	return alignment, dbs


def parse_seed_alignments(ali_directory, dbs_directory):
	alignments = list()
	dbss = list()
	filenames = os.listdir(ali_directory)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists(dbs_directory + filename_base + '.dbs'):
			alignment, dbs = parse_seed_alignment(ali_directory, dbs_directory, filename_base)
			alignments.append(alignment)
			dbss.append(dbs)
		else:
			print('Warning: No dbs data found for \'' + filename + '\'.')
	return alignments, dbss


def parse_alignments_dataset(ali_set_path, dbs_path, filename):
	with open(ali_set_path + filename + '.alis') as file:
		alignments = list()
		ali_idx = -1
		line = file.readline()
		while line != '':
			content = line.split()

			if content[1].isnumeric():
				ali_idx += 1
				alignments.append(list())
			else:
				alignments[ali_idx].append(content[1])

			line = file.readline()

	with open(dbs_path + filename + '.dbs') as dbs_file:
		dbs = dbs_file.readlines()[1].split()[0]

	return alignments, dbs


def parse_alignments_datasets(ali_sets_directory, dbs_directory):
	alignments_sets = list()
	dbss = list()
	filenames = os.listdir(ali_sets_directory)
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists(dbs_directory + filename_base + '.dbs'):
			alignments, dbs = parse_alignments_dataset(ali_sets_directory, dbs_directory, filename_base)
			alignments_sets.append(alignments)
			dbss.append(dbs)
		else:
			print('Warning: No dbs data found for \'' + filename + '\'.')

	alignments_sets_consecutively = list()
	dbss_consecutively = list()
	for alignments_set, dbs in zip(alignments_sets, dbss):
		for alignment in alignments_set:
			alignments_sets_consecutively.append(alignment)
			dbss_consecutively.append(dbs)
	return alignments_sets_consecutively, dbss_consecutively


def main():
	seed_alignments, seed_dbss = parse_seed_alignments('../rnaconv/data/Rfam.seed_frequency/ali/',
									  '../rnaconv/data/Rfam.seed_neighbourhood/dbs/')
	sissi_alignment_datasets, sissi_dbss = parse_alignments_datasets('../rnaconv/data/sissi/',
														 '../rnaconv/data/Rfam.seed_neighbourhood/dbs/')


if __name__ == '__main__':
	main()
