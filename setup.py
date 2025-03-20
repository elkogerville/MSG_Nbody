from setuptools import setup, find_packages

setup(
    name='MSG_Nbody',
    version='0.1.1',
    description='Nbody simulation code for galaxy interactions',
    url='https://github.com/elkogerville/MSG_Nbody',
    author='Elko Gerville-Reache',
    author_email='elkogerville@gmail.com',
    license='MIT',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'numpy',
        'numba',
        'tqdm',
        'matplotlib',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
