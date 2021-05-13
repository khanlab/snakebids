import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snakebids",
    version="0.3.10",
    author="Ali Khan",
    author_email="alik@robarts.ca",
    description="BIDS integration into snakemake workflows",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/akhanf/snakebids",
    include_package_data=True,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "snakebids-create=snakebids.template:main",
        ],
    },
    install_requires=[
        "pybids==0.12.3",
        "snakemake>=5.28.0",
        "PyYAML>=5.3.1",
        "cookiecutter>=1.7.2",
    ],
    python_requires=">=3.7",
)
