from setuptools import setup, find_packages

VERSION = '0.5'

long_description = ''

setup(
    name='gridgen',
    version=VERSION,
    description='',
    long_description=long_description,
    url='https://github.com/gpiantoni/gridgen',
    author="Gio Piantoni",
    author_email='gridgen@gpiantoni.com',
    license='BSD',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
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
            'gridgen=gridgen.bin.command:main',
        ],
    },
    )
