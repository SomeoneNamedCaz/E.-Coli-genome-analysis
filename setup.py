from setuptools import *

packageDir = "src"


setup(
	name="snpalign",
	version="0.0.10",
	description="a pseudo-multiple alignment program and end to end GWAS tool",
	package_dir={"":packageDir},
	packages=find_packages(where=packageDir, exclude=['tests', ".git"]),
	entry_points={'console_scripts': [
        'snpalign = snpalign.endToEndScript:main',
    ]},
	python_requires=">=3.7",
)
