import os


def parse_alignments(filepath, filename):
	alignments = list()
	ali_idx = -1
	with open(filepath + filename) as file:
		line = file.readline()
		while line != '':
			content = line.split()

			if content[1].isnumeric():
				ali_idx += 1
				alignments.append(list())
			else:
				alignments[ali_idx].append(content[1])

			line = file.readline()

	return alignments


def parse_alignment_files(directory):
	filenames = os.listdir(directory)
	data = list()
	for filename in filenames:
		alignments = parse_alignments(directory, filename)
		data.append(alignments)
	return data


def main():
	data = parse_alignment_files('../rnaconv/data/sissi/')
	pass


if __name__ == '__main__':
	main()
