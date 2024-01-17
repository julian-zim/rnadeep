import os


def main():
	ali_path = './rfam/seed_alignments/'
	freq_path = './rfam/seed_frequencies/'
	tree_path = './rfam/seed_trees/fixed/'

	filenames = os.listdir(ali_path)
	for filename in filenames:
		keep = False
		with open(ali_path + filename) as file:
			line = file.readline()
			while line != '':
				for char in line:
					if char in ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']:
						keep = True
						break
				line = file.readline()
		if not keep:
			os.remove(ali_path + filename)
			del_single_freq_path = freq_path + 'single/' + filename.split('.')[0] + '.freq'
			if os.path.exists(del_single_freq_path):
				os.remove(del_single_freq_path)
			del_double_freq_path = freq_path + 'doublet/' + filename.split('.')[0] + '.freq'
			if os.path.exists(del_double_freq_path):
				os.remove(del_double_freq_path)
			del_tree_path = tree_path + filename.split('.')[0] + '.seed_tree'
			if os.path.exists(del_tree_path):
				os.remove(del_tree_path)


if __name__ == '__main__':
	main()
