from setuptools import setup

setup(name = 'ptmonte',
      version = '0.1',
      description = 'Monte Carlo simulation package',
      url = 'https://github.com/zeldery/ptmonte',
      author = 'Thien-Phuc Tu-Nguyen',
      licence = 'GNU',
      packages = ['ptmonte'],
      install_requires = ['numpy','pandas'],
      include_package_data = True,
      zip_safe = False)
