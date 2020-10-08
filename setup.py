import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()

setuptools.setup(
    name='rip_multilevel',
    version='0.1',
    scripts=['rip_multilevel'],
    authos="Evangelos A. Dimopoulos, Evan K. Irving-Pease",
    author_email='antonisdim41@gmail.com',
    description='Species dientification pipeline for both single species and metagenomic samples.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/antonisdim/rip",
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 3",
                 "License :: MIT License"],
)


