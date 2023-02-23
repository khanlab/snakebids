import json

import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("pipeline_description.json", "r", encoding="utf-8") as fh:
    pipeline = json.load(fh)
name = pipeline["GeneratedBy"][0]["Name"]
description = pipeline["Name"]
version = pipeline["GeneratedBy"][0]["Version"]
optional_vals = {
    attr: pipeline["GeneratedBy"][0].get(key)
    for attr, key in [
        ("url", "CodeURL"),
        ("author", "Author"),
        ("author_email", "AuthorEmail"),
    ]
    if key in pipeline
}

setuptools.setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    long_description_content_type="text/x-rst",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "{{cookiecutter.app_name}}={{cookiecutter.app_name}}.run:main"
        ]
    },
    install_requires=["snakebids>={{cookiecutter.snakebids_version}}", "snakemake"],
    python_requires=">=3.7",
    **optional_vals,
)
