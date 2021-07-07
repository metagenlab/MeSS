from setuptools import setup, find_packages
from mess import __version__

setup(name='MeSS',
      version=__version__,
      url='https://github.com/metagenlab/MeSS',
      license='GPL-3',
      author='Farid Chaabane',
      author_email='faridchaabn@gmail.com',
      description="MeSS is a snakemake workflow used for simulating metagenomic mock communities",
      packages=find_packages(),
      include_package_data=True,
      data_files=[(".", ["LICENSE", "README.md"])],
      install_requires=['click', 'snakemake==6.2.1'],
      entry_points={
            'console_scripts': [
                  'mess = mess.mess:cli'
            ]
      },
      )