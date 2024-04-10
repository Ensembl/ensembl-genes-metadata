# pylint: disable=missing-module-docstring
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import argparse
import re
import json
import os
import ast
from typing import Dict, Union

# Add options for each key
keys_to_extract = {
    "started_job": "Started job on",
    "started_mapping": "Started mapping on",
    "finished": "Finished on",
    "mapping_speed": "Mapping speed, Million of reads per hour",
    "number_of_input_reads": "Number of input reads",
    "average_input_read_length": "Average input read length",
    "uniquely_mapped_reads_number": "Uniquely mapped reads number",
    "uniquely_mapped_reads_percentage": "Uniquely mapped reads %",
    "average_mapped_length": "Average mapped length",
    "number_of_splices_total": "Number of splices: Total",
    "number_of_splices_annotated": "Number of splices: Annotated (sjdb)",
    "number_of_splices_GT_AG": "Number of splices: GT/AG",
    "number_of_splices_GC_AG": "Number of splices: GC/AG",
    "number_of_splices_AT_AC": "Number of splices: AT/AC",
    "number_of_splices_non_canonical": "Number of splices: Non-canonical",
    "mismatch_rate_per_base": "Mismatch rate per base, %",
    "deletion_rate_per_base": "Deletion rate per base",
    "deletion_average_length": "Deletion average length",
    "insertion_rate_per_base": "Insertion rate per base",
    "insertion_average_length": "Insertion average length",
    "number_of_reads_mapped_to_multiple_loci": "Number of reads mapped to multiple loci",
    "percentage_reads_mapped_to_multiple_loci": "% of reads mapped to multiple loci",
    "number_of_reads_mapped_to_too_many_loci": "Number of reads mapped to too many loci",
    "percentage_reads_mapped_to_too_many_loci": "% of reads mapped to too many loci",
    "number_of_reads_unmapped_too_many_mismatches": "Number of reads unmapped: too many mismatches",
    "percentage_reads_unmapped_too_many_mismatches": "% of reads unmapped: too many mismatches",
    "number_of_reads_unmapped_too_short": "Number of reads unmapped: too short",
    "percentage_reads_unmapped_too_short": "% of reads unmapped: too short",
    "number_of_reads_unmapped_other": "Number of reads unmapped: other",
    "percentage_reads_unmapped_other": "% of reads unmapped: other",
    "number_of_chimeric_reads": "Number of chimeric reads",
    "percentage_of_chimeric_reads": "% of chimeric reads",
}


def parse_star_output(file_path:str, keys: dict, output_dir:str, extra_parameters: Dict[str, str])-> str:
    """Parse STAR Log file

    Args:
        file_path (str): path Log.out
        keys (dict): keys to inspect
        output_dir (str): Path outout dir
        extra_parameters (Dict[str, str]): Extra info to add to the json

    Returns:
        str: _description_
    """

    result = {}
    for parameter, value in extra_parameters.items():
        result[parameter] = value
    with open(file_path, "r") as file:
        for line in file:
            for key in keys:
                match = re.search(f"{key}\s+\|\s+(.+)", line) #pylint: disable=anomalous-backslash-in-string
                if match:
                    result[key] = match.group(1).strip()
    output_file = os.path.join(output_dir, "star_output.json")
    with open(output_file, "w") as json_file:
        json.dump(result, json_file, indent=4)
    return output_file


def parse_extra_parameters(param_str: Union[str, ast.AST]):
    """Check type Dict for extra parameters

    Args:
        param_str Union[str, ast.AST]): extra_parameters 

    Raises:
        argparse.ArgumentTypeError: exception wrong format

    Returns:
        _type_: type
    """
    try:
        return ast.literal_eval(param_str)
    except ValueError:
        raise argparse.ArgumentTypeError(#pylint:disable=raise-missing-from
            "Invalid extra parameters format. Must be a valid Python dictionary."
        )


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse STAR output and create JSON with specified keys")
    parser.add_argument("--file_path", required=True, type=str, help="Path to the STAR output file")
    parser.add_argument("--output_dir", required=True, type=str, help="Output dir")
    parser.add_argument(
        "--extra_parameters", required=False, type=parse_extra_parameters, help="{'key':'value'}"
    )
    for key in keys_to_extract:
        parser.add_argument(f"--{key}", action="store_true", help=f"Include {key} in output JSON")
    # print(vars(parser))
    return parser.parse_args()


def main():
    """Parser entry-point."""
    args = parse_args()

    options = vars(args)
    # del options['file_path']
    # del options['run_accession']
    # del options['output_dir']
    # if the corresponding key is present in the options dictionary
    # (meaning the user has provided the corresponding command-line argument)
    keys_to_include = [value for key, value in keys_to_extract.items() if options.get(key)]

    output_json = parse_star_output(args.file_path, keys_to_include, args.output_dir, args.extra_parameters)
    return output_json


if __name__ == "__main__":
    main()
