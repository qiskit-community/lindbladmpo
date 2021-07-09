# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import os
import re
import platform
import subprocess

from setuptools import setup, find_packages, Extension

with open('requirements') as f:
	REQUIREMENTS = f.read().splitlines()

s_system = platform.system().lower()
if s_system == 'windows':
	s_executable = 'lindbladmpo.exe'
	print('Note: Automatic building of the solver executable is not currently supported on Windows.\n'
		  '\tPlease follow the installation guide in the pacakge directory for instructions regarding\n'
		  '\tlocall installation of cygwin, and building the solver from source.')
else:
	s_executable = 'lindbladmpo.exe'
	process = subprocess.Popen('make -C ./src/', shell=True)
	exit_code = process.wait()
	process = subprocess.Popen(f'cp ./bin/{s_executable} ./lindbladmpo/{s_executable}', shell=True)

# Read long description from README.
README_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)),
						   'README.md')
with open(README_PATH) as readme_file:
	README = re.sub(
		'<!--- long-description-skip-begin -->.*<!--- long-description-skip-end -->', '',
		readme_file.read(), flags=re.S | re.M)

setup(
	name="lindblad-mpo",
	version="0.1.0",
	description="A matrix-product-operators solver for the Lindblad master equation",
	long_description=README,
	long_description_content_type='text/markdown',
	url="https://github.com/haggaila/lindbladmpo",
	author="Gregoire Misguich",
	author_email="gregoire.misguich@ipht.fr",
	license="Apache 2.0",
	classifiers=[
		"Environment :: Console",
		"License :: OSI Approved :: Apache Software License",
		"Intended Audience :: Developers",
		"Intended Audience :: Science/Research",
		# "Operating System :: Microsoft :: Windows",
		"Operating System :: MacOS",
		"Operating System :: POSIX :: Linux",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Topic :: Scientific/Engineering",
	],
	keywords="mps mpo lindblad sdk quantum",
	packages=find_packages(exclude=['test*']),
	install_requires=REQUIREMENTS,
	include_package_data=True,
	python_requires=">=3.6",
	extras_require={
		'visualization': ['matplotlib>=2.1'],
	},
	project_urls={
		"Source Code": "https://github.com/haggaila/lindbladmpo",
	},
	zip_safe=False
)
