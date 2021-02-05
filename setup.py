from setuptools import setup
from mess import __version__

setup(name='MeSS',
      version=__version__,
      url='https://github.com/metagenlab/MeSS',
      license='BSD-3',
      author='Farid Chaabane',
      description="MeSS is a snakemake workflow used for simulating metagenomic mock communities",
      packages=['mess'],
      scripts=['mess/scripts/add_read_percent.py',
               'mess/scripts/krona_table.py',
               'mess/scripts/merge_contigs.py',
               'mess/scripts/remove_index.py',
               'mess/scripts/shuffle.sh',
               'mess/scripts/simulate_reads.rules',
               'mess/scripts/Snakefile',
               'mess/assembly_finder/Snakefile',
               'mess/assembly_finder/rules/assembly_table.py',
               'mess/assembly_finder/rules/combine_tables.py',
               'mess/assembly_finder/rules/dl.py',
               'mess/assembly_finder/rules/filter_table.py',
               'mess/assembly_finder/rules/find_assemblies.rules'],
      package_data={'': [
            "mess/*",
      ]},
      data_files=[(".", ["README.md"])],
      include_package_data=True,
      install_requires=['click==7', 'snakemake-minimal==5.32'],
      entry_points={
            'console_scripts': [
                  'mess = mess.mess:cli'
            ]
      },
      )