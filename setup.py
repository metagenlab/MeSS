import os
from setuptools import setup, find_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "mess",
            "mess.VERSION",
        )
    ) as f:
        return f.readline().strip()


def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="mess",
    packages=find_packages(),
    url="https://github.com/metagenlab/MeSS",
    python_requires=">=3.11",
    description="Snakemake pipeline for simulating short and long read metagenomes",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Farid Chaabane",
    author_email="farid.chaabane@chuv.ch",
    data_files=get_data_files(),
    py_modules=["mess"],
    install_requires=[
        "snakemake>=8.0.0,<=8.24.1",
        "snaketool-utils>=0.0.5",
        "attrmap>=0.0.7",
        "pandas>=2.2.1",
        "biopython>=1.83",
        "rich-click>=1.8.3",
    ],
    entry_points={"console_scripts": ["mess=mess.__main__:main"]},
    include_package_data=True,
)
