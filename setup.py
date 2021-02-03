from setuptools import setup


def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]

setup(name='MeSS',
      version=get_version("mess/__init__.py"),
      url='https://github.com/metagenlab/MeSS',
      license='BSD-3',
      author='Farid Chaabane',
      description="MeSS is a snakemake workflow used for simulating metagenomic mock communities",
      )