from setuptools import setup, find_packages

setup(
    name='MSG_Nbody',
    version='0.2.5',
    description='Nbody simulation code for galaxy interactions',
    url='https://github.com/elkogerville/MSG_Nbody',
    author='Elko Gerville-Reache',
    author_email='elkogerville@gmail.com',
    license='MIT',
    long_description=open('README_MSG.md').read(),
    long_description_content_type='text/markdown',
    include_package_data=True,
    packages=find_packages(exclude=[
        'ANIMATIONS*',
        'DOCUMENTATION*',
        'Initial_Conditions*',
        'Tests*',
	'MSG_Nbody.egg-info*',
        '*.egg-info*',
        'build*',
        'dist*',
    ]),
    install_requires=[
        'numpy',
        'numba',
        'tqdm',
        'matplotlib',
        'scipy',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
