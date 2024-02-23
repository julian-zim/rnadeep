import sys
import os
import shutil


def copy(from_path, to_path, filenames):
	tree_orig_path = 'seed_trees/original'
	tree_fixed_path = 'seed_trees/fixed'
	tree_rescaled_path = 'seed_trees/rescaled'

	ali_path = 'seed_alignments'

	neigh_wuss_path = 'seed_neighbourhoods/wuss'
	neigh_dbn_path = 'seed_neighbourhoods/dbn'
	neigh_ct_path = 'seed_neighbourhoods/ct'
	neigh_nei_path = 'seed_neighbourhoods/nei'

	single_freq_path = 'seed_frequencies/single'
	doublet_freq_path = 'seed_frequencies/doublet'

	os.makedirs(os.path.join(to_path, tree_orig_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, tree_fixed_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, tree_rescaled_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, ali_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_wuss_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_dbn_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_ct_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, neigh_nei_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, single_freq_path), exist_ok=True)
	os.makedirs(os.path.join(to_path, doublet_freq_path), exist_ok=True)

	for filename in filenames:
		condition = [os.path.exists(os.path.join(from_path, tree_orig_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, tree_rescaled_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, ali_path, filename + '.aln')),
				os.path.exists(os.path.join(from_path, neigh_wuss_path, filename + '.wuss')),
				os.path.exists(os.path.join(from_path, neigh_dbn_path, filename + '.dbn')),
				os.path.exists(os.path.join(from_path, neigh_ct_path, filename + '.ct')),
				os.path.exists(os.path.join(from_path, neigh_nei_path, filename + '.nei')),
				os.path.exists(os.path.join(from_path, single_freq_path, filename + '.freq')),
				os.path.exists(os.path.join(from_path, doublet_freq_path, filename + '.freq'))]
		if all(condition):
			shutil.copy(os.path.join(from_path, tree_orig_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_orig_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_fixed_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, tree_rescaled_path, filename + '.seed_tree'),
						os.path.join(to_path, tree_rescaled_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(from_path, ali_path, filename + '.aln'),
						os.path.join(to_path, ali_path, filename + '.aln'))
			shutil.copy(os.path.join(from_path, neigh_wuss_path, filename + '.wuss'),
						os.path.join(to_path, neigh_wuss_path, filename + '.wuss'))
			shutil.copy(os.path.join(from_path, neigh_dbn_path, filename + '.dbn'),
						os.path.join(to_path, neigh_dbn_path, filename + '.dbn'))
			shutil.copy(os.path.join(from_path, neigh_ct_path, filename + '.ct'),
						os.path.join(to_path, neigh_ct_path, filename + '.ct'))
			shutil.copy(os.path.join(from_path, neigh_nei_path, filename + '.nei'),
						os.path.join(to_path, neigh_nei_path, filename + '.nei'))
			shutil.copy(os.path.join(from_path, single_freq_path, filename + '.freq'),
						os.path.join(to_path, single_freq_path, filename + '.freq'))
			shutil.copy(os.path.join(from_path, doublet_freq_path, filename + '.freq'),
						os.path.join(to_path, doublet_freq_path, filename + '.freq'))
		else:
			print('Skipping \'' + filename + '\' as it is missing data. (Raw output: ' + str(condition) + ')')


def main():
	if len(sys.argv) < 4:
		print('Usage: ./copy.py <path-to-copy-from> <path-to-copy-to> <filenames-to-copy> [<additional-filenames-to-copy> ...]')
		return -1

	copy(sys.argv[1], sys.argv[2], sys.argv[3:])


if __name__ == '__main__':
	main()
