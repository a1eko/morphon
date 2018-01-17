from setuptools import setup

setup(name='morphon',
      version='0.0.1',
      description='neuron morphology analysis tool',
      author='Alexander Kozlov',
      author_email='akozlov@kth.se',
      license='BSD',
      packages=['morphon'],
      scripts=[
          'apps/mview',
          'apps/mfeatures',
          'apps/mrepair',
      ],
      install_requires=[
          'numpy',
      ],
      zip_safe=False)
