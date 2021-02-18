import setuptools

from haystac import __version__

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt") as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

setuptools.setup(
    name="haystac",
    version=__version__,
    author="Evangelos A. Dimopoulos, Evan K. Irving-Pease",
    author_email="antonisdim41@gmail.com",
    description="Species identification pipeline for both single species and metagenomic samples.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/antonisdim/haystac",
    package_dir={
        "haystac": "haystac",
        "haystac.workflow": "haystac/workflow",
        "haystac.workflow.scripts": "haystac/workflow/scripts",
    },
    packages=["haystac", "haystac.workflow", "haystac.workflow.scripts"],
    include_package_data=True,
    entry_points={"console_scripts": ["haystac = haystac.cli:Haystac"]},
    classifiers=["Programming Language :: Python :: 3", "License :: OSI Approved :: MIT License"],
    license="MIT",
    install_requires=required,
)
