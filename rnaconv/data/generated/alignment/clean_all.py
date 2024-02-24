import os
import subprocess


def clean(full=True):
	datasets = list(set(os.listdir())
					- {'__pycache__'}
					- set([subdir for subdir in os.listdir() if '.' in subdir.split('/')[-1]])
					- {'full'} if full else set())

	for dataset in datasets:
		command = 'cd ./' + dataset + '/; ./clean.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def main():
	clean(full=False)


if __name__ == '__main__':
	main()
