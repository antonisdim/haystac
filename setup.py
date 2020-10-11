import setuptools
from glob import glob

from haystack.workflow import __version__, _program

with open('README.md', 'r') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

setuptools.setup(
    name='haystack_metagenomics',
    version='0.1',
    author="Evangelos A. Dimopoulos, Evan K. Irving-Pease",
    author_email='antonisdim41@gmail.com',
    description='Species identification pipeline for both single species and metagenomic samples.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/antonisdim/rip",
    package_dir = {'haystack': 'haystack',
            'haystack.workflow': 'haystack/workflow',
            'haystack.workflow.scripts': 'haystack/workflow/scripts'},
    packages=['haystack', 'haystack.workflow', 'haystack.workflow.scripts'],
    include_package_data=True,
    entry_points="""
      [console_scripts]
      {program} = haystack.haystack:Rip
      """.format(program = _program),
    classifiers=["Programming Language :: Python :: 3",
                 "License :: MIT License"],
    license='MIT',
    install_requires=required,
)
