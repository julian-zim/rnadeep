import os
import pdfkit

# sudo apt-get install wkhtmltopdf; sudo apt-get install python-dev-is-python3; pip install pdfkit

directories = ['../../examples', '../../rnaconv', '../../rnadeep']
html_files = []
for directory in directories:
	for filename in os.listdir(directory):
		if filename.endswith('.py'):
			filepath = os.path.join(directory, filename)
			os.system(f'pydoc -w {filepath}')
			html_file = os.path.basename(filepath).replace('.py', '.html')
			html_files.append(html_file)

merged_html = ''
for html_file in html_files:
	try:
		with open(html_file, 'r') as f:
			merged_html += f.read()
	except FileNotFoundError:
		print('Skipping ' + html_file + ', as it\'s not existing.')

pdfkit.from_string(merged_html, '../documentation.pdf')
