from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    'biopython>=1.70',
]

desc = """A collection of functions for subsetting and reformatting sequences from fasta files."""

setup(name='fastsub',
      version='1.0.0',
      description=desc,
      url='https://github.com/Adamtaranto/subset-fasta-by-name',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['fastsub'],
      classifiers=pypi_classifiers,
      keywords=["fasta"],
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'fastsub=fastsub.cmd_line:main',
        ],
    },
    )
