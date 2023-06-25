from __future__ import annotations

from shutil import rmtree
from sys import version_info as py_version

if py_version >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata


def update_files(files: list[str] | str, replacement_str: str, cc_variable: str):
    """Helper function to update cookiecutter content for certain files

    INPUT
    ------
    files
        file(s) to be updated
    replacement_str
        text to substitute into file(s)
    cc_var
        cookiecutter variable to replace (e.g.
        {{ cookiecutter._snakebids_version }})
    """
    for fpath in files:
        with open(fpath, "r", encoding="utf-8") as fcontent:
            content = fcontent.read()

        content = content.replace(cc_variable, replacement_str)

        with open(fpath, "w", encoding="utf-8") as fcontent:
            fcontent.write(content)


# Replace snakebids version in cookiecutter
sb_file_lists = ["setup.py", "docs/requirements.txt"]
sb_version = metadata.version("snakebids")
update_files(sb_file_lists, sb_version, "{{ cookiecutter._snakebids_version }}")

# Remove documentation template
if not {{cookiecutter.create_doc_template}}:  # noqa: F821
    rmtree("docs")
