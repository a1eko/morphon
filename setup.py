import setuptools


def readme():
    with open('README.rst') as f:
        return f.read()


setuptools.setup(name="morphon",
    version="0.0.4",
    description="Neuron morphology processing scripts",
    long_description=readme(),
    author="Alexander Kozlov",
    author_email="akozlov@kth.se",
    url="https://github.com/a1eko/morphon",
    packages=setuptools.find_packages(),
    scripts=[
        "apps/mcheck",
        "apps/mfind",
        "apps/mmeter",
        "apps/mmod",
        "apps/mrep",
        "apps/msholl",
        "apps/mview",
    ],
    install_requires=[
        "numpy",
        "matplotlib",
    ],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX",
    ],
    keywords=(
        "neuron",
        "morphology",
        "reconstruction",
        "morphometry",
        "modification",
        "repair"),
)
