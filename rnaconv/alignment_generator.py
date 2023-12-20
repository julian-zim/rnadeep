def wuss_to_ct(secondary_structure):
	pass


# stockholm format
def get_neighbourshoods(alignment):
	pass


# stockholm format
def get_frequencies(alignment):
	pass


# stockholm format
def convert_stockholm_alignment(filename):
	return '', '', ''


def get_parameters(neighbourhood_file, frequencies_file, tree_file):
	return [], [], []


def generate_alignments(parameters, n):
	pass


def main():
	n = 10

	filename = ''
	filenames = convert_stockholm_alignment(filename)
	parameters = get_parameters(*filenames)
	generate_alignments(parameters, n)


if __name__ == '__main__':
	main()
