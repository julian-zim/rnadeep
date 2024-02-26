import sys
import os


def filter_rfam_data(ali_dirpath,
					 single_freq_dirpath,
					 doublet_freq_dirpath,
					 neigh_wuss_dirpath,
					 neigh_dbn_dirpath,
					 tree_fixed_dirpath,
					 tree_rescaled_dirpath,
					 max_length=700):
	"""
	Filters the alignments, consensus structures, frequencies and trees of a converted rfam database (not touching
	the original tree files and original Rfam.seed file):
	Every data point (consisting of the four filetypes mentioned above) which alignment exceeds a certain length is
	removed.

		Parameters:
			ali_dirpath (str): Path to the directory of the converted alignments.
			single_freq_dirpath (str): Path to the directory of the extracted single frequency files.
			doublet_freq_dirpath (str): Path to the directory of the extracted doublet frequency files.
			neigh_wuss_dirpath (str): Path to the directory of the extracted wuss files.
			neigh_dbn_dirpath (str): Path to the directory of the extracted dbn files.
			tree_fixed_dirpath (str): Path to the directory of the fixed newick string tree files.
			tree_rescaled_dirpath (str): Path to the directory of the rescaled tree files.
			max_length (int, None, optional): Maximum allowed length of an alignment. Default is 700.
	"""

	if max_length is None:
		max_length = 700

	for filename in os.listdir(ali_dirpath):
		ali_filepath = os.path.join(ali_dirpath, filename)
		with open(ali_filepath) as file:
			file.readline()
			length = len(file.readline().split()[1])

		filename_base = filename.split('.')[0]
		if length > max_length:
			os.remove(ali_filepath)
			try:
				os.remove(os.path.join(single_freq_dirpath, filename_base + '.sfreq'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(doublet_freq_dirpath, filename_base + '.dfreq'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(neigh_wuss_dirpath, filename_base + '.wuss'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(neigh_dbn_dirpath, filename_base + '.dbn'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(tree_fixed_dirpath, filename_base + '.seed_tree'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(tree_rescaled_dirpath, filename_base + '.seed_tree'))
			except FileNotFoundError:
				pass

			print('Removed data of alignment \'' + filename_base + '\' due to it exceeding a length of ' + str(max_length) + '.')


def main():
	if not 2 <= len(sys.argv) <= 3:
		print('Usage: ./rfam_filter.py <rfam_directory_path> [max_length]')
		return -1

	rfam_dirpath = sys.argv[1]
	max_length = None
	if len(sys.argv) == 3:
		max_length = int(sys.argv[2])

	ali_filepath = os.path.join(rfam_dirpath, 'seed_alignments')
	single_freq_filepath = os.path.join(rfam_dirpath, 'seed_frequencies/single')
	doublet_freq_filepath = os.path.join(rfam_dirpath, 'seed_frequencies/doublet')
	neigh_wuss_filepath = os.path.join(rfam_dirpath, 'seed_neighbourhoods/wuss')
	neigh_dbn_filepath = os.path.join(rfam_dirpath, 'seed_neighbourhoods/dbn')
	treefixed_filepath = os.path.join(rfam_dirpath, 'seed_trees/fixed')
	treerescaled_filepath = os.path.join(rfam_dirpath, 'seed_trees/rescaled')
	if (not os.path.exists(ali_filepath)
			or not os.path.exists(single_freq_filepath)
			or not os.path.exists(doublet_freq_filepath)
			or not os.path.exists(neigh_wuss_filepath)
			or not os.path.exists(neigh_dbn_filepath)
			or not os.path.exists(treefixed_filepath)
			or not os.path.exists(treerescaled_filepath)):
		print('Warning: The given path doesn\'t look like a converted Rfam database. We will filter what we find.\n'
			  'Expected folder structure:\n'
			  '- rfam\n\t- seed_alignments\n\t- seed_frequencies\n\t\t- doublet\n\t\t- single\n\t- seed_neighbourhoods'
			  '\n\t\t- dbn\n\t\t- wuss\n\t- seed_trees\n\t\t- fixed\n\t\t- original\n\t\t- rescaled')

	print('========== FILTERING RFAM DATA POINTS ==========')
	filter_rfam_data(ali_filepath,
					 single_freq_filepath,
					 doublet_freq_filepath,
					 neigh_wuss_filepath,
					 neigh_dbn_filepath,
					 treefixed_filepath,
					 treerescaled_filepath,
					 max_length)
	print('Done.')


if __name__ == '__main__':
	main()
