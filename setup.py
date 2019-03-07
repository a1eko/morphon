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
          'apps/mfind',
          'apps/mmeter',
          'apps/mmod',
          'apps/mrep',
          'apps/msholl',
          'apps/mshow',
      ],
      install_requires=[
          'numpy',
          'matplotlib',
      ],
      zip_safe=False)
