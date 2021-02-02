from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='BO_interfaces',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/SaeedM75/BO_interfaces',
    packages=['BO_interfaces'],
    install_requires=['numpy', 'ase==3.19.1', 'pymatgen', 'bayesian-optimization', 'pandas', 'scikit-learn==0.21.3'],
)
