from setuptools import setup, find_packages

with open('README.md', 'r') as readme:
    long_description = readme.read()

setup(
    name='pyWavelet-TBC',
    version='0.0.1',
    author='Climate Research Division',
    author_email='',
    description='A package for working with ',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='path_to_github_repo',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: canada_open?',
        'Operating System :: OS Independent',
    ],
    setup_requires=['pytest-runner', ],
    tests_requires=['pytest', ],
)