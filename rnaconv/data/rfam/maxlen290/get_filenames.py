import os


def main():
	ali_path = '../full/seed_alignments'

	filenames = os.listdir(os.path.join(ali_path))
	ambig_filenames = list()
	for filename_full in filenames:
		with open(os.path.join(ali_path, filename_full)) as file:
			file.readline()
			if len(file.readline().split()[1]) <= 290:
				ambig_filenames.append(filename_full.split('.')[0])

	print(' '.join(ambig_filenames))


if __name__ == '__main__':
	main()
