import sys
import os


default_max_length = 700


def filter_rfam_data(rfam_path, max_length):
	if max_length is None:
		global default_max_length
		max_length = default_max_length

	treerescaled_filepath = os.path.join(rfam_path, 'seed_trees/fixed')
	treefixed_filepath = os.path.join(rfam_path, 'seed_trees/fixed')
	ali_filepath = os.path.join(rfam_path, 'seed_alignments')

	neigh_nei_filepath = os.path.join(rfam_path, 'seed_neighbourhoods/nei')
	neigh_ct_filepath = os.path.join(rfam_path, 'seed_neighbourhoods/ct')
	neigh_dbn_filepath = os.path.join(rfam_path, 'seed_neighbourhoods/dbn')
	neigh_wuss_filepath = os.path.join(rfam_path, 'seed_neighbourhoods/wuss')

	single_freq_filepath = os.path.join(rfam_path, 'seed_frequencies/single')
	doublet_freq_filepath = os.path.join(rfam_path, 'seed_frequencies/doublet')

	for filename in os.listdir(ali_filepath):
		with open(os.path.join(ali_filepath, filename)) as file:
			file.readline()
			length = len(file.readline().split()[1])

		filename_base = filename.split('.')[0]
		if length > max_length:
			os.remove(os.path.join(ali_filepath, filename))

			try:
				os.remove(os.path.join(treefixed_filepath, filename_base + '.seed_tree'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(treerescaled_filepath, filename_base + '.seed_tree'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(neigh_nei_filepath, filename_base + '.nei'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(neigh_ct_filepath, filename_base + '.ct'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(neigh_dbn_filepath, filename_base + '.dbn'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(neigh_wuss_filepath, filename_base + '.wuss'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(single_freq_filepath, filename_base + '.sfreq'))
			except FileNotFoundError:
				pass
			try:
				os.remove(os.path.join(doublet_freq_filepath, filename_base + '.dfreq'))
			except FileNotFoundError:
				pass

			print('Removed data of alignment \'' + filename_base + '\' due to it exceeding a length of ' + str(max_length) + '.')


def main():
	if len(sys.argv) < 2:
		print('Usage: ./rfam_converter.py <rfam_path> [max_length]')
		return -1

	rfam_path = sys.argv[1]

	max_length = None
	if len(sys.argv) == 3:
		max_length = int(sys.argv[2])

	print('========== FILTERING RFAM DATA POINTS ==========')
	filter_rfam_data(rfam_path, max_length)
	print('Done')


if __name__ == '__main__':
	main()
