import os
import subprocess


def generate_alignment(filename, n):
	base_path = './data/'
	#ct_filepath = base_path + 'Rfam.seed_neighbourhood/ct/' + filename + 'fa.ct'
	neigh_filepath = base_path + 'Rfam.seed_neighbourhood/nei/' + filename + '.fa.ct.nei'
	freq_filepath = base_path + 'Rfam.seed_frequency/freq/' + filename + '.fa.freq'
	tree_filepath = base_path + 'Rfam.seed_tree/fixed/' + filename + '.seed_tree'
	out_filepath = base_path + 'sissi/' + filename + '_' + str(n) + '.ali'

	with open(neigh_filepath) as nei_file:
		seq_length = sum(1 for _ in nei_file)
	with open(freq_filepath) as freq_file:
		frequencies = freq_file.readline()[:-1]
	command = ('./sissi099' +
			  #' -nr' + ct_filepath +
			   ' -nn ' + neigh_filepath +
			   ' -fs ' + frequencies +
			   ' -l' + str(seq_length) +
			   ' -a' + str(n) +
			   ' ' + tree_filepath +
			   ' > ' + out_filepath)
	result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if os.path.getsize(out_filepath) == 0:
		#print('Running \'' + command + '\'')
		print('Error in \'' + filename + '\':\n' + result.stderr)
		os.remove(out_filepath)
	else:
		print('Successfully generated alignment(s) for \'' + filename + '\'.')


def main():
	try:
		os.makedirs('./data/sissi/')
	except FileExistsError:
		pass

	filenames = os.listdir('./data/Rfam.seed_neighbourhood/nei/')
	for filename in filenames:
		filename_base = filename.split('.')[0]
		if os.path.exists('./data/Rfam.seed_tree/fixed/' + filename_base + '.seed_tree'):
			generate_alignment(filename_base, 10)


if __name__ == '__main__':
	main()
