from setuptools import setup

setup(name='morphon',
      version='0.0.2',
      description='neuron morphology processing tools',
      author='Alexander Kozlov',
      author_email='akozlov@kth.se',
      license='GPLv3',
      packages=['morphon'],
      scripts=[
          'apps/mcheck',
          'apps/mshow',
          'apps/mfind',
          'apps/mrep',
          'apps/mmod',
          'apps/mmeter',
      ],
      install_requires=[
          'numpy',
          'matplotlib',
      ],
      zip_safe=False)
