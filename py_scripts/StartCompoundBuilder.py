import argparse
import json
import logging
import os
import sys
import time
from typing import List

import requests

import BuildCompoundDir
"""
Warning! this script will fail unless it has its requirements.txt requirements installed!
"""

class ReadableDir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir = values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace, self.dest, prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


def main(args):
    cal_default_dest = '/Users/cmartin/projects/fake_compound_dir/'
    cal_ftp = '/Users/cmartin/projects/fakeftppub/'
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w', '--destination', action=ReadableDir, help="Output directory", default=cal_default_dest)
    parser.add_argument('-f', '--ftp', action=ReadableDir, default=cal_ftp, help="FTP directory")
    parser.add_argument('-n', '--new_compounds_only', action="store_true", help="whether to only run the builder for new chebi entries")
    args = parser.parse_args(args)

    destination_directory = args.destination
    ftp = args.ftp
    new_compounds_only = args.new_compounds_only
    overall_start_time = time.time()

    with open(f'{ftp}reactome.json') as reac_file:
        reactome_data = json.load(reac_file)

    with open(f'{ftp}mapping.json') as map_file:
        ml_mapping = json.load(map_file)

    compound_list = []
    try:
        compound_list = requests.get("http://www.ebi.ac.uk/metabolights/webservice/compounds/list").json()['content']
    except Exception as e:
        logging.exception('unable to get compounds list.')
        exit(0)

    if new_compounds_only:
        dircomps = get_mtblc_ids_from_directory(destination_directory)
        delta = get_delta(webservice_list=compound_list, filesystem_list=dircomps)
        print(delta)
        compound_list = delta

    for compound in compound_list:
        current_compound_start_time = time.time()
        # build process returns dict, no use for it in prod but handy when debugging
        __ = BuildCompoundDir.build(
            metabolights_id=compound.strip(), ml_mapping=ml_mapping,reactome_data=reactome_data,
            data_directory=destination_directory)
        current_compound_end_time = time.time()
        print(f'time taken to process {compound}: {current_compound_end_time - current_compound_start_time}')

    overall_end_time = time.time()
    diff_delta = time.gmtime(overall_end_time - overall_start_time)
    hours, minutes, seconds = diff_delta.tm_hour, diff_delta.tm_min, diff_delta.tm_sec
    print(f'Time taken for compound building process: {hours:02d}:{minutes:02d}:{seconds:02d}')


def get_mtblc_ids_from_directory(directory: str):
    list_of_ids = []
    for entry in os.scandir(directory):
        if entry.is_dir():
            list_of_ids.append(entry.name)
    return list_of_ids


def get_delta(webservice_list: List[str], filesystem_list: List[str]) -> List[str]:
    return list(set(webservice_list) - set(filesystem_list))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
