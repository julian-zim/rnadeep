import sys
import os
import shutil


def copy(out_path,
		 ali_dirpath,
		 single_freq_dirpath,
		 doublet_freq_dirpath,
		 neigh_wuss_dirpath,
		 neigh_dbn_dirpath,
		 tree_fixed_dirpath,
		 tree_rescaled_dirpath,
		 filenames):
	"""
	Copy selected alignments, consensus structures, frequencies and trees of a converted rfam database to a new location.

		Parameters:
			out_path (str): Path to which to copy the data to
			ali_dirpath (str): Path to the directory of the converted alignment files
			single_freq_dirpath (str): Path to the directory of the extracted single frequency files
			doublet_freq_dirpath (str): Path to the directory of the extracted doublet frequency files
			neigh_wuss_dirpath (str): Path to the directory of the extracted wuss files
			neigh_dbn_dirpath (str): Path to the directory of the extracted dbn files
			tree_fixed_dirpath (str): Path to the directory of the fixed newick string tree files
			tree_rescaled_dirpath (str): Path to the directory of the rescaled tree files
			filenames (list): List of filenames (without ending) to copy the alignments, frequencies, secondary
				structures and trees of
	"""

	ali_outpath = os.path.join(out_path, 'seed_alignments')
	single_freq_outpath = os.path.join(out_path, 'seed_frequencies/single')
	doublet_freq_outpath = os.path.join(out_path, 'seed_frequencies/doublet')
	neigh_wuss_outpath = os.path.join(out_path, 'seed_neighbourhoods/wuss')
	neigh_dbn_outpath = os.path.join(out_path, 'seed_neighbourhoods/dbn')
	tree_fixed_outpath = os.path.join(out_path, 'seed_trees/fixed')
	tree_rescaled_outpath = os.path.join(out_path, 'seed_trees/rescaled')

	os.makedirs(ali_outpath, exist_ok=True)
	os.makedirs(single_freq_outpath, exist_ok=True)
	os.makedirs(doublet_freq_outpath, exist_ok=True)
	os.makedirs(neigh_wuss_outpath, exist_ok=True)
	os.makedirs(neigh_dbn_outpath, exist_ok=True)
	os.makedirs(tree_fixed_outpath, exist_ok=True)
	os.makedirs(tree_rescaled_outpath, exist_ok=True)

	for filename in filenames:
		condition = [os.path.exists(os.path.join(ali_dirpath, filename + '.aln')),
					 os.path.exists(os.path.join(single_freq_dirpath, filename + '.sfreq')),
					 os.path.exists(os.path.join(doublet_freq_dirpath, filename + '.dfreq')),
					 os.path.exists(os.path.join(neigh_wuss_dirpath, filename + '.wuss')),
					 os.path.exists(os.path.join(neigh_dbn_dirpath, filename + '.dbn')),
					 os.path.exists(os.path.join(tree_fixed_dirpath, filename + '.seed_tree')),
					 os.path.exists(os.path.join(tree_rescaled_dirpath, filename + '.seed_tree'))]
		if all(condition):
			shutil.copy(os.path.join(ali_dirpath, filename + '.aln'),
						os.path.join(ali_outpath, filename + '.aln'))
			shutil.copy(os.path.join(single_freq_dirpath, filename + '.sfreq'),
						os.path.join(single_freq_outpath, filename + '.sfreq'))
			shutil.copy(os.path.join(doublet_freq_dirpath, filename + '.dfreq'),
						os.path.join(doublet_freq_outpath, filename + '.dfreq'))
			shutil.copy(os.path.join(neigh_wuss_dirpath, filename + '.wuss'),
						os.path.join(neigh_wuss_outpath, filename + '.wuss'))
			shutil.copy(os.path.join(neigh_dbn_dirpath, filename + '.dbn'),
						os.path.join(neigh_dbn_outpath, filename + '.dbn'))
			shutil.copy(os.path.join(tree_fixed_dirpath, filename + '.seed_tree'),
						os.path.join(tree_fixed_outpath, filename + '.seed_tree'))
			shutil.copy(os.path.join(tree_rescaled_dirpath, filename + '.seed_tree'),
						os.path.join(tree_rescaled_outpath, filename + '.seed_tree'))
		else:
			print('Skipping \'' + filename + '\' as it is missing data. (Raw output: ' + str(condition) + ')')


def main():
	if len(sys.argv) < 4:
		print('Usage: ./copy.py <path_to_copy_to> <rfam_path_to_copy_from> <filenames_to_copy> '
			  '[<additional_filenames_to_copy> ...]')
		return -1

	outpath = sys.argv[1]
	rfam_path = sys.argv[2]
	filenames = sys.argv[3:]

	ali_dirpath = os.path.join(rfam_path, 'seed_alignments')
	single_freq_dirpath = os.path.join(rfam_path, 'seed_frequencies/single')
	doublet_freq_dirpath = os.path.join(rfam_path, 'seed_frequencies/doublet')
	neigh_wuss_dirpath = os.path.join(rfam_path, 'seed_neighbourhoods/wuss')
	neigh_dbn_dirpath = os.path.join(rfam_path, 'seed_neighbourhoods/dbn')
	tree_fixed_dirpath = os.path.join(rfam_path, 'seed_trees/fixed')
	tree_rescaled_dirpath = os.path.join(rfam_path, 'seed_trees/rescaled')
	if (not os.path.exists(ali_dirpath)
			or not os.path.exists(single_freq_dirpath)
			or not os.path.exists(doublet_freq_dirpath)
			or not os.path.exists(neigh_wuss_dirpath)
			or not os.path.exists(neigh_dbn_dirpath)
			or not os.path.exists(tree_fixed_dirpath)
			or not os.path.exists(tree_rescaled_dirpath)):
		print('Warning: The given path doesn\'t look like a converted Rfam database. We will copy what we find.\n'
			  'Expected folder structure:\n'
			  '- rfam\n\t- seed_alignments\n\t- seed_frequencies\n\t\t- doublet\n\t\t- single\n\t- seed_neighbourhoods'
			  '\n\t\t- dbn\n\t\t- wuss\n\t- seed_trees\n\t\t- fixed\n\t\t- rescaled')

	copy(outpath,
		 ali_dirpath,
		 single_freq_dirpath,
		 doublet_freq_dirpath,
		 neigh_wuss_dirpath,
		 neigh_dbn_dirpath,
		 tree_fixed_dirpath,
		 tree_rescaled_dirpath,
		 filenames)
	print('Done.')


if __name__ == '__main__':
	main()
