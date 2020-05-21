from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
   name='itex_codes',
   version='1.0',
   description='Codes to create and mantain CII db',
   license='GNU',
   long_description=long_description,
   author='Eric March Vila',
   author_email='eric.march@upf.edu',
   url='https://github.com/phi-grib/Itex_codes',
   packages=['CreateDB','UpdateDB']
)
