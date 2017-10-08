from setuptools import setup

setup(
    name='medi_gwas',
    version='0.0.1',
    description=
    'Maps relationships between biomedical data, \
        compares genetic drug-phenotype pairs (GWAS) with clinical database (MEDI)',
    long_description='',
    packages=['medi_gwas'],
    install_requires=[],
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'pandacli=panda.panda:cli',
            ]},
    author='Wen Yin',
    author_email='wen.yin@duke.edu',
    classifiers=(
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5'
        ),
)


if __name__ == '__main__':
    pass