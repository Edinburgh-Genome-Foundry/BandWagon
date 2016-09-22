import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('gelsimulator/version.py').read()) # loads __version__

setup(name='gelsimulator',
      version=__version__,
      author='Zulko',
      description='',
      long_description=open('README.rst').read(),
      license='see LICENSE.txt',
      keywords="Gel Agarose Simulation Matplotlib",
      packages= find_packages(exclude='docs')
)
