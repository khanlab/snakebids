import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snakebids", 
    version="0.1.1",
    author="Ali Khan",
    author_email="alik@robarts.ca",
    description="BIDS integration into snakemake workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/akhanf/snakebids",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pybids>=0.12.0",
        "snakemake>=5.25",
        "PyYAML>=5.3",
    ],
    python_requires='>=3.7'
)
