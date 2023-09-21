#!/usr/bin/env python


import datetime
import sys
import os
import json
import re
import logging
import argparse

import libchebipy

from libchebipy._chebi_entity import ChebiEntity


class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


def configure_logger() -> logging.Logger:
    logger = logging.getLogger('valLogger')

    # create handlers - one for stdout, one to file
    #s_handler = logging.StreamHandler()
    f_handler = logging.FileHandler(f'val_{datetime.datetime.now()}.log')
   # s_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.INFO)

    # set formatters for each handler
   # s_format = logging.Formatter('%(name)s // %(levelname)s // %(message)s')
    f_format = logging.Formatter('%(asctime)s // %(funcName)s // %(name)s // %(levelname)s // %(message)s')
    #s_handler.setFormatter(s_format)
    f_handler.setFormatter(f_format)

    #logger.addHandler(s_handler)
    logger.addHandler(f_handler)
    return logger


def main(args):
    # all this script does is go over the results of MLCompoundsBot and report any discrepancies. It does not alter
    # the reference directories. It is primarily to check whether MTBLC12345 correponds to CHEBI12345, but also checks if
    # the MTBLC files are in the right format, and counts how many non MTBLC files there are, and how many entities libchebipy
    # was unable to find.
    cal_default_dest = '/Users/cmartin/projects/fake_compound_dir/'
    ebi_default_dest = '/nfs/www-prod/web_hx2/cm/metabolights/prod/reference/'
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', '--compound_dir', action=readable_dir,help="directory where the compound subdirs live", default=cal_default_dest)
    args = parser.parse_args(args)
    dir_path = args.compound_dir

    logger = configure_logger()
    highest = len([name for name in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, name))])


    # counters for where we don't care about tracking individual IDs
    correct_counter = 0
    non_mtblc_counter = 0

    # registers for when we want to inspect individual MTBLC objects after
    incorrect_register = []
    bad_egg_register = []
    chebi_no_found_register = []
    just_bad_ids = []

    # walk through every subdirectory in the reference directory
    for root, dirs, files in os.walk(dir_path, topdown=False):
        # skip if subdir is empty
        if len(files) == 0:
            logger.warning(f'No Files found in {root}')
            non_mtblc_counter += 1
            continue
        # skip if not a root MTBLC subdirectory
        if files[0].startswith('MTBLC') is False:
            non_mtblc_counter += 1
            continue
        with open(f'{root}/{files[0]}') as mtblc_json:
            data = json.load(mtblc_json)
            # if we don't have an ID its not an MTBLC json
            if 'id' not in data.keys():
                logger.warning(f'{files[0]} is not valid mtblc json')
                bad_egg_register.append(files[0])
                continue
            ostensible_chebi_id = re.split('([0-9]+)', data['id'])[1]
            try:
                chebi_entity = ChebiEntity(ostensible_chebi_id)
            except libchebipy.ChebiException as e:
                logging.exception(f'Couldnt retrieve chebi entity: {str(e)}')
                chebi_no_found_register.append(files[0])

            chebi_name = chebi_entity.get_name()
            # check if the chebi and mtblc names match - may need to do some fuzzy matching here.
            if chebi_name.lower() == data['name'].lower():
                logger.warning(f'{files[0]} is correct')
                correct_counter += 1
            else:
                logger.warning(f'{files[0]} is incorrect')
                incorrect_register.append(f'{files[0]} / chebi name:{chebi_name} [{data["id"]}] / mtblc name: {data["name"]} [{ostensible_chebi_id}]\n')
                just_bad_ids.append(f'{ostensible_chebi_id}\n')
    # output the offending ones to a text file here
    with open(dir_path + 'incorrect.txt', 'w') as incorrect:
        incorrect.writelines(incorrect_register)
    with open(dir_path + 'bad_ids.txt', 'w') as bad_ids:
        bad_ids.writelines(just_bad_ids)
    logger.warning(f'difference: {highest - correct_counter} \n '
                f'bad files: {len(bad_egg_register)} \n '
                f'non mtblc files: {non_mtblc_counter} \n'
                f'chebi couldnt find: {len(chebi_no_found_register)}')


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
