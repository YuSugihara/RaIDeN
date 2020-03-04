#!/usr/bin/env python

from distutils.core import setup
from raiden.__init__ import __version__

setup(name='raiden',
      version='{}'.format(__version__),
      description='RaIDeN: Rapid identification of causal gene with target motif using NGS technology',
      author='Yu Sugihara',
      author_email='yu57th@gmail.com',
      url='https://github.com/YuSugihara/RaIDeN',
      license='GPL',
      packages=['raiden'],
      entry_points={'console_scripts': [
            'raiden = raiden.raiden:main',
            'jiji = raiden.jiji:main',
            ]
        }
     )