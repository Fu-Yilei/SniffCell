from setuptools import setup, find_packages
from src.__init__ import __version__
setup(
    name='sniffcell',
    version=__version__,
    packages=find_packages(),
    url='https://github.com/Fu-Yilei/SniffCell',
    license='MIT',
    author='Yilei Fu',
    author_email='yilei.fu@bcm.edu',
    description='SniffCell: Annotate SVs cell type based on CpG methylation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)