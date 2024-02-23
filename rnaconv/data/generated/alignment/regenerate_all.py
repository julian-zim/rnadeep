import os
import subprocess


def generate(preclean=True, full=True):
	datasets = list(set(os.listdir())
					- {'__pycache__'}
					- set([subdir for subdir in os.listdir() if '.' in subdir.split('/')[-1]])
					- {'full'})
	if full:
		datasets.insert(0, 'full')

	for dataset in datasets:
		if preclean:
			command = 'cd ./' + dataset + '/; ./clean.sh'
			print('Running \'' + command + '\'')
			subprocess.run(command, shell=True, text=True)

		command = 'cd ./' + dataset + '/; ./generate.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def filtering(full=True):
	datasets = list(set(os.listdir())
					- {'__pycache__'}
					- set([subdir for subdir in os.listdir() if '.' in subdir.split('/')[-1]])
					- {'full'})
	if full:
		datasets.insert(0, 'full')

	for dataset in datasets:
		command = 'cd ./' + dataset + '/; ./filter.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def main():
	generate(preclean=True, full=True)
	filtering(full=True)


if __name__ == '__main__':
	main()
