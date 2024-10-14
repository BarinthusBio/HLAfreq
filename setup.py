from setuptools import setup, find_packages

setup(
    name="HLAfreq",
    version="0.0.5",
    url="https://github.com/BarinthusBio/HLAfreq",
    project_urls={
        'Documentaion': "https://barinthusbio.github.io/HLAfreq/HLAfreq.html",
        'Tracker': "https://github.com/BarinthusBio/HLAfreq/issues"
    },
    author="David Wells",
    author_email="david.wells@BarinthusBio.co.uk",
    description="Download and combine HLA frequency data from multiple studies",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    python_requires='>=3.10',
    install_requires=[
        'bs4',
        'requests',
        'pandas',
        'numpy>=1.24.0',
        'matplotlib>=3.6.0',
        'scipy>=1.10.0',
        'pymc>=3',
        'arviz'
    ],
    classifiers=[
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
    ],
    package_dir={'':'src'},
    packages=find_packages(
        where='src'
    ),
)