import os
import subprocess


def convert(preclean=True, full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		if preclean:
			command = 'cd ./rfam/full/; ./clean.sh'
			print('Running \'' + command + '\'')
			subprocess.run(command, shell=True, text=True)

		command = 'cd ./rfam/full/; ./convert.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)

	for dataset in datasets:
		if preclean:
			command = 'cd ./rfam/' + dataset + '/; ./clean.sh'
			print('Running \'' + command + '\'')
			subprocess.run(command, shell=True, text=True)

		command = 'cd ./rfam/' + dataset + '/; ./copy.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def rfam_filter(full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		command = 'cd ./rfam/full/; ./filter.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)

	for dataset in datasets:
		command = 'cd ./rfam/' + dataset + '/; ./filter.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def generate(preclean=True, full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		if preclean:
			command = 'cd ./sissi/full/; ./clean.sh'
			print('Running \'' + command + '\'')
			subprocess.run(command, shell=True, text=True)

		command = 'cd ./sissi/full/; ./generate.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)

	for dataset in datasets:
		if preclean:
			command = 'cd ./sissi/' + dataset + '/; ./clean.sh'
			print('Running \'' + command + '\'')
			subprocess.run(command, shell=True, text=True)

		command = 'cd ./sissi/' + dataset + '/; ./generate.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def sissi_filter(full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		command = 'cd ./sissi/full/; ./filter.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)

	for dataset in datasets:
		command = 'cd ./sissi/' + dataset + '/; ./filter.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def main():
	convert(preclean=True, full=True)
	rfam_filter(full=True)

	generate(preclean=True, full=True)
	sissi_filter(full=True)


if __name__ == '__main__':
	main()
