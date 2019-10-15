import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open("bandwagon/version.py").read())  # loads __version__

setup(
    name="bandwagon",
    version=__version__,
    author="Zulko",
    description="Simulate DNA band patterns for gel migration experiments",
    url="https://github.com/Edinburgh-Genome-Foundry/BandWagon",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="Gel Agarose Simulation Matplotlib",
    packages=find_packages(exclude="docs"),
    install_requires=(
        "biopython",
        "matplotlib",
        "numpy",
        "scipy",
        "snapgene_reader",
    ),
)
