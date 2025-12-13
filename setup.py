from setuptools import setup, find_packages

setup(
    name='sniffcell',
    version='0.5.3',
    packages=find_packages(),
    url='https://github.com/Fu-Yilei/SniffCell',
    license='MIT',
    author='Yilei Fu',
    author_email='yilei.fu@bcm.edu',
    description='SniffCell: Annotate SVs cell type based on CpG methylation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)