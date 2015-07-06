#!/usr/bin/env python

from distutils.core import setup

# Get version number from package
exec(open('readgeneratorlib/version.py').read())

setup(
    name='ReadGenerator',
    version=__version__,
    description='Creates synthetic paired-end DNA reads from given fraction of Human, Bacteria, Phix174 and Virus/Phage.',
    author='Ashwini Patil',
    author_email='patil.ashwini1091@gmail.com',
    url='https://github.com/apatil1/',
    packages=['readgeneratorlib'],
    scripts=['readgenerator.py'],
    )
