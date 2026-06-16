from setuptools import setup, find_packages

setup(
    name='bio_kit',
    version='1.0.0',
    description='Bioinformatics analysis toolkit for teaching and small-scale research',
    author='Bio Kit Team',
    packages=find_packages(exclude=['web', 'tests', 'web.*', 'tests.*']),
    include_package_data=True,
    install_requires=[
        'biopython>=1.78',
        'numpy>=1.20.0',
        'matplotlib>=3.4.0',
        'ete3>=3.1.2',
    ],
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
    ],
)
