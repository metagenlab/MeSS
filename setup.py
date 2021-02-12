from setuptools import setup,find_packages
from mess import __version__

setup(name='MeSS',
      version=__version__,
      url='https://github.com/metagenlab/MeSS',
      license='GPL-3',
      author='Farid Chaabane',
      author_email='faridchaabn@gmail.com',
      description="MeSS is a snakemake workflow used for simulating metagenomic mock communities",
      packages=find_packages(),
      scripts=['mess/scripts/populate_table_percent.py',
               'mess/scripts/read_counts_table.py',
               'mess/scripts/merge_contigs.py',
               'mess/scripts/shuffle.sh',
               'mess/scripts/simulate_reads.rules',
               'mess/scripts/Snakefile',
               'mess/assembly_finder/Snakefile',
               'mess/assembly_finder/rules/assembly_table.py',
               'mess/assembly_finder/rules/combine_tables.py',
               'mess/assembly_finder/rules/dl.py',
               'mess/assembly_finder/rules/filter_table.py',
               'mess/assembly_finder/rules/find_assemblies.rules',
               "mess/assembly_finder/envs/Assembly_finder.yml",
               "mess/envs/simulate_reads.yml"],
      include_package_data=True,
      data_files=[(".", ["LICENSE.txt","README.md"])],
      install_requires=['click', 'snakemake'],
      entry_points={
            'console_scripts': [
                  'mess = mess.mess:cli'
            ]
      },
      )