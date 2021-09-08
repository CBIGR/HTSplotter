from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
  name = 'HTSplotter',         # How you named your package folder (MyLib)
  packages = ['HTSplotter'],   # Chose the same as "name"
  version = '0.14',      # Start with a small number and increase it with every change you make
  license='GNU General Public License v3 or later (GPLv3+)',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'An end-to-end data processing, '
                'analysis and visualisation tool for chemical and genetic in '
                'vitro perturbation screens',   # Give a short description about your library
  long_description=long_description,
  long_description_content_type="text/markdown",
  author = 'CarolinadCNunes',                   # Type in your name
  author_email = 'carolina.decarvalhonunes@ugent.be',      # Type in your E-Mail
  url = 'https://github.com/CBIGR/HTSplotter',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/CBIGR/HTSplotter/archive/refs/tags/HTSplotter_v0.14.tar.gz',    # I explain this later on
  keywords = ['High-throughput screening', 'drug combination',
              'genetic-chemical perturbation', 'dose-response'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
                    'certifi',
                    'cycler',
                    'h5py',
                    'kiwisolver',
                    'matplotlib',
                    'minio',
                    'numpy',
                    'Pillow',
                    'psutil',
                    'pyparsing',
                    'pyPdf',
                    'PyPDF2',
                    'PyPDF3',
                    'python-dateutil',
                    'pytz',
                    'scipy',
                    'seaborn',
                    'six',
                    'tqdm',
                    'urllib3',
                    'xlrd',
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
  ],
)
