import os
import subprocess


def convert(preclean=True, full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		if preclean:
			command = 'cd ./rfam/full/; ./clean.sh'
			print('Running \'' + command + '\'')
			process = subprocess.Popen(command, shell=True, text=True)
			process.wait()

		command = 'cd ./rfam/full/; ./convert.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()

	for dataset in datasets:
		if preclean:
			command = 'cd ./rfam/' + dataset + '/; ./clean.sh'
			print('Running \'' + command + '\'')
			process = subprocess.Popen(command, shell=True, text=True)
			process.wait()

		command = 'cd ./rfam/' + dataset + '/; ./copy.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()


def rfam_filter(full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		command = 'cd ./rfam/full/; ./filter.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()

	for dataset in datasets:
		command = 'cd ./rfam/' + dataset + '/; ./filter.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()


def generate(preclean=True, full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		if preclean:
			command = 'cd ./sissi/full/; ./clean.sh'
			print('Running \'' + command + '\'')
			process = subprocess.Popen(command, shell=True, text=True)
			process.wait()

		command = 'cd ./sissi/full/; ./generate.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()

	for dataset in datasets:
		if preclean:
			command = 'cd ./sissi/' + dataset + '/; ./clean.sh'
			print('Running \'' + command + '\'')
			process = subprocess.Popen(command, shell=True, text=True)
			process.wait()

		command = 'cd ./sissi/' + dataset + '/; ./generate.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()


def sissi_filter(full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		command = 'cd ./sissi/full/; ./filter.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()

	for dataset in datasets:
		command = 'cd ./sissi/' + dataset + '/; ./filter.sh'
		print('Running \'' + command + '\'')
		process = subprocess.Popen(command, shell=True, text=True)
		process.wait()


def regenerate(preclean=True, rfam_filtering=True, sissi_filtering=True, full=True):
	convert(preclean, full)
	if rfam_filtering:
		rfam_filter(full)

	generate(preclean, full)
	if sissi_filtering:
		sissi_filter(full)


def main():
	regenerate(True, True, True, True)


if __name__ == '__main__':
	main()
