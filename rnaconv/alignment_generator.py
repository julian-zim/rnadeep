import os
import subprocess


def generate_alignments(outpath):
	try:
		os.makedirs(outpath)
	except FileExistsError:
		pass

	filenames = os.listdir('./data/Rfam.seed_neighbourhood/nei/')
	for filename in filenames:
		filename_base = filename.split('.')[0]
		with open('./data/Rfam.seed_frequency/freq/' + filename_base + '.fa.freq') as freq_file:
			frequencies = freq_file.readline()[:-1]
		with open('./data/Rfam.seed_neighbourhood/nei/' + filename_base + '.fa.ct.nei') as nei_file:
			seq_length = sum(1 for _ in nei_file)
		command = ('./sissi099' +
				  #' -nr./data/Rfam.seed_neighbourhood/ct/' + filename_base + '.fa.ct' +
				   ' -nn ./data/Rfam.seed_neighbourhood/nei/' + filename_base + '.fa.ct.nei' +
				   ' -fs ' + frequencies +
				   ' -l' + str(seq_length) +
				   ' -a1' +
				   ' ./data/Rfam.seed_tree/fixed/' + filename_base + '.seed_tree' +
				   ' > ./data/sissi/' + filename_base + '.ali')
		result = subprocess.run(command, text=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		if os.path.getsize('./data/sissi/' + filename_base + '.ali') == 0:
			# print('Running \'' + command + '\'')
			print('Error for file \'' + filename_base + '.ali\':\n' + result.stderr)
			os.remove('./data/sissi/' + filename_base + '.ali')


def main():
	generate_alignments('./data/sissi/')


if __name__ == '__main__':
	main()
