from setuptools import setup, find_packages


setup(
    name='cheminformatics',
    version='0.0.3',
    license='MIT',
    author='Jacob Gerlach',
    author_email='jwgerlach00@gmail.com',
    url='https://github.com/jwgerlach00/cheminformatics',
    description='Data science and machine learning tools for chemistry applications',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    python_requires='>=3.6',
    install_requires=[
        'naclo',
        'numpy',
        'rdkit',
        'rdkit_pypi',
        'scikit_learn'
    ],
)
