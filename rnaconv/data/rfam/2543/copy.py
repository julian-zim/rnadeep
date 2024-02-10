import os
import shutil


def copy(filenames):
	to_path = '.'
	from_path = '../full'

	tree_orig_path = 'seed_trees/original'
	tree_fixed_path = 'seed_trees/fixed'

	ali_path = 'seed_alignments'

	neigh_wuss_path = 'seed_neighbourhoods/wuss'
	neigh_dbn_path = 'seed_neighbourhoods/dbn'
	neigh_ct_path = 'seed_neighbourhoods/ct'
	neigh_nei_path = 'seed_neighbourhoods/nei'

	single_freq_path = 'seed_frequencies/single'
	doublet_freq_path = 'seed_frequencies/doublet'

	try:
		os.makedirs(os.path.join(os.path.join(to_path, tree_orig_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, tree_fixed_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, ali_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, neigh_wuss_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, neigh_dbn_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, neigh_ct_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, neigh_nei_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, single_freq_path)))
	except FileExistsError:
		pass
	try:
		os.makedirs(os.path.join(os.path.join(to_path, doublet_freq_path)))
	except FileExistsError:
		pass

	for filename in filenames:
		condition = [os.path.exists(os.path.join(from_path, tree_orig_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree')),
				os.path.exists(os.path.join(from_path, ali_path, filename + '.ali')),
				os.path.exists(os.path.join(from_path, neigh_wuss_path, filename + '.wuss')),
				os.path.exists(os.path.join(from_path, neigh_dbn_path, filename + '.dbn')),
				os.path.exists(os.path.join(from_path, neigh_ct_path, filename + '.ct')),
				os.path.exists(os.path.join(from_path, neigh_nei_path, filename + '.nei')),
				os.path.exists(os.path.join(from_path, single_freq_path, filename + '.freq')),
				os.path.exists(os.path.join(from_path, doublet_freq_path, filename + '.freq'))]
		if all(condition):
			shutil.copy(os.path.join(os.path.join(from_path, tree_orig_path, filename + '.seed_tree')),
									 os.path.join(to_path, tree_orig_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(os.path.join(from_path, tree_fixed_path, filename + '.seed_tree')),
									 os.path.join(to_path, tree_fixed_path, filename + '.seed_tree'))
			shutil.copy(os.path.join(os.path.join(from_path, ali_path, filename + '.ali')),
									 os.path.join(to_path, ali_path, filename + '.ali'))
			shutil.copy(os.path.join(os.path.join(from_path, neigh_wuss_path, filename + '.wuss')),
									 os.path.join(to_path, neigh_wuss_path, filename + '.wuss'))
			shutil.copy(os.path.join(os.path.join(from_path, neigh_dbn_path, filename + '.dbn')),
									 os.path.join(to_path, neigh_dbn_path, filename + '.dbn'))
			shutil.copy(os.path.join(os.path.join(from_path, neigh_ct_path, filename + '.ct')),
									 os.path.join(to_path, neigh_ct_path, filename + '.ct'))
			shutil.copy(os.path.join(os.path.join(from_path, neigh_nei_path, filename + '.nei')),
									 os.path.join(to_path, neigh_nei_path, filename + '.nei'))
			shutil.copy(os.path.join(os.path.join(from_path, single_freq_path, filename + '.freq')),
									 os.path.join(to_path, single_freq_path, filename + '.freq'))
			shutil.copy(os.path.join(os.path.join(from_path, doublet_freq_path, filename + '.freq')),
									 os.path.join(to_path, doublet_freq_path, filename + '.freq'))


def main():
	filenames = ['RF02543']
	copy(filenames)


if __name__ == '__main__':
	main()
