import json
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('pipeline_description.json', 'r') as fh:
    pipeline = json.load(fh)
    print(pipeline)
    name = pipeline['GeneratedBy'][0]['Name']
    description = pipeline['Name']
    version = pipeline['GeneratedBy'][0]['Version']
    url = pipeline['GeneratedBy'][0]['CodeURL']
    

setuptools.setup(
    name=name,
    version=version,
    author="Ali Khan",
    author_email="alik@robarts.ca",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=url,
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pybids==0.12.3",
        "snakemake==5.28.0",
        "PyYAML==5.3.1",
    ],
    python_requires='>=3.7'
)
