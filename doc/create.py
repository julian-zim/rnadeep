import os
import pdfkit

# sudo apt-get install wkhtmltopdf; sudo apt-get install python-dev-is-python3; pip install pdfkit

directories = ['../../examples', '../../rnaconv', '../../rnadeep']
excluded_files = ['lstm_models.html', 'metrics.html', 'mlforensics.html', 'sliding_window.html']

html_files = []
for directory in directories:
	for filename in os.listdir(directory):
		if filename.endswith('.py'):
			filepath = os.path.join(directory, filename)
			os.system(f'pydoc -w {filepath}')
			html_file = os.path.basename(filepath).replace('.py', '.html')
			if html_file not in excluded_files:
				html_files.append(html_file)

merged_html = ''
for html_file in html_files:
	try:
		with open(html_file, 'r') as f:
			merged_html += f.read()
	except FileNotFoundError:
		print('Skipping ' + html_file + ', as it\'s not existing.')

filtered_html = merged_html.replace('/home/julian-zim/Files/Cloud/OneDrive/OneFiles/Linux/Work/Workspaces/Study/UNIVIE/PyCharm/PR_SPB/','')

pdfkit.from_string(filtered_html, '../documentation.pdf')
