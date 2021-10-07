""" Script to fix dataset_description.json PyBIDS compatibility """
from pathlib import Path
import json
import argparse


def get_parser():
    parser = argparse.ArgumentParser(
        description="Update BIDS Derivatives Dataset  to be backwards compatible with pybids",
        epilog=("Adds PipelineDescription.Name to the dataset_description.json file, from earlier BIDS spec) that pybids requires, using the first entry of the GeneratedBy list.  This is needed if parsing a modern version of fmriprep derivatives with snakebids/pybids."),
    )
    parser.add_argument("dataset_root", help="BIDS derivative dataset to fix")
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    root_path = Path(args.dataset_root)

    # check if a dataset_description.json file exists
    json_path = root_path.joinpath("dataset_description.json")

    # read the dataset_description json file
    with open(json_path) as f:
        dataset_dict = json.load(f)

    if "PipelineDescription" in dataset_dict.keys():
        if "Name" in dataset_dict["PipelineDescription"].keys():
            print(f"{json_path} already compatible with pybids")
            return

    pipeline_name = dataset_dict["GeneratedBy"][0]["Name"]
    pipeline_version = dataset_dict["GeneratedBy"][0]["Version"]

    # add to the dict
    print(
        f"Adding PipelineDescription.Name: {pipeline_name}, PipelineDescription.Version: {pipeline_version}"
    )

    dataset_dict["PipelineDescription"] = {
        "Name": pipeline_name,
        "Version": pipeline_version,
    }

    print(f"Saving to {json_path}")
    with open(json_path, "w") as write_file:
        json.dump(dataset_dict, write_file, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()
