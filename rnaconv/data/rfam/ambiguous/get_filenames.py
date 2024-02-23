import os


def main():
	ali_path = '../full/seed_alignments'

	filenames = os.listdir(os.path.join(ali_path))
	ambig_filenames = list()
	for filename_full in filenames:
		add = False
		with open(os.path.join(ali_path, filename_full)) as file:
			line = file.readline()
			while line != '':
				split_line = line.split()
				if len(split_line) > 1:
					for char in split_line[1]:
						if char in ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']:
							add = True
							break
				line = file.readline()
		if add:
			ambig_filenames.append(filename_full.split('.')[0])

	print(' '.join(ambig_filenames))


if __name__ == '__main__':
	main()
