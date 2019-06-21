from setuptools import setup
setup(
	name='pysamiterators',
	version='0.01',
	author='Buys de Barbanson',
	author_email='code@buysdb.nl',
	description='Pysam related iterators',
	url='https://github.com/BuysDB/pysamiterators',
	packages=['pysamiterators',
        'pysamiterators.iterators'
        ],
	install_requires=['pysam']
)
