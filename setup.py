from setuptools import setup
setup(
	name='pysamiterators',
	version='1.08',
	author='Buys de Barbanson',
	author_email='code@buysdb.nl',
	description='Pysam related iterators',
	url='https://github.com/BuysDB/pysamiterators',
	download_url='https://github.com/BuysDB/pysamiterators/archive/v1.08.tar.gz',
	packages=['pysamiterators','pysamiterators.iterators'
        ],
	install_requires=['pysam', 'numpy'],
	 keywords = ['mate-pair','pysam','iteration'],
	license='MIT',
	classifiers=[
	 'Development Status :: 4 - Beta',
	  'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
	  'License :: OSI Approved :: MIT License'
	]
)
