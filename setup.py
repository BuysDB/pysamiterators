from setuptools import setup
setup(
	name='pysamiterators',
	version='1.0',
	author='Buys de Barbanson',
	author_email='code@buysdb.nl',
	description='Pysam related iterators',
	url='https://github.com/BuysDB/pysamiterators',
	download_url='https://github.com/BuysDB/pysamiterators/archive/v1.0.tar.gz',
	packages=['pysamiterators'
        ],
	install_requires=['pysam'],
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
