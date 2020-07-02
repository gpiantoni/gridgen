from setuptools import setup, find_packages

VERSION = '0.4'

long_description = ''

setup(
    name='gridloc',
    version=VERSION,
    description='',
    long_description=long_description,
    url='https://github.com/gpiantoni/gridloc',
    author="Gio Piantoni",
    author_email='gridloc@gpiantoni.com',
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.6',
        ],
    keywords='analysis',
    packages=find_packages(exclude=('test', )),
    install_requires=[
        'numpy',
        'scipy',
        'wonambi',
        'nibabel',
        'plotly',
        ],
    entry_points={
        'console_scripts': [
            'gridloc=gridloc.bin.command:main',
        ],
    },
    )
