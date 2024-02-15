import os
import subprocess


def clean_rfam(full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		command = 'cd ./rfam/full/; ./clean.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)

	for dataset in datasets:
		command = 'cd ./rfam/' + dataset + '/; ./clean.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def clean_sissi(full=True):
	datasets = list(set(os.listdir('rfam')) - {'full'})

	if full:
		command = 'cd ./sissi/full/; ./clean.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)

	for dataset in datasets:
		command = 'cd ./sissi/' + dataset + '/; ./clean.sh'
		print('Running \'' + command + '\'')
		subprocess.run(command, shell=True, text=True)


def clean(full=True):
	clean_rfam(full)
	clean_sissi(full)


def main():
	clean(True)


if __name__ == '__main__':
	main()
