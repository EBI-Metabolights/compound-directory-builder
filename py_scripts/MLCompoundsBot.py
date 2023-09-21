#!/usr/bin/env python

# Author: Venkata chandrasekhar Nainala 
# Version: 0.1.0
# Email: mailcs76@gmail.com / venkata@ebi.ac.uk
# Date: 19 May 2016

""" 
    Dependencies:
        Python 2.7
"""
import concurrent.futures
import sys
import argparse
import threading

import requests

import utils
import logging
import os
import time
import json
import subprocess
from random import randint

destinationDirectory = ""
workingDirectory = ""
mlSCMappingFile = ""
reactomeJSONFile = ""
ftp = ""
importLevel = 0

reactomeData = {}
mlMapping = {}
requestedCompound = ""

ebi_default_dest = '/nfs/www-prod/web_hx2/cm/metabolights/prod/reference/'
cal_default_dest = '/Users/cmartin/projects/fake_compound_dir/'

ebi_ftp = '/net/isilonP/public/rw/homes/tc_cm01/from_ftp_pub/compounds/'
cal_ftp = '/Users/cmartin/projects/fakeftppub/'

compound_thread = threading.local()
file_lock = threading.Lock()

def get_session():
    if not hasattr(compound_thread, 'session'):
        compound_thread.session = requests.Session()
    return compound_thread.session

def arg_wrapper(compound):
    session = get_session()
    logging.info("-----------------------------------------------")
    logging.info("Fetching compound: " + compound)
    return utils.fetchCompound(
        compound.strip(), workingDirectory, destinationDirectory, reactomeData, mlMapping, session, file_lock)

def get_all_compounds(compoundList):
    # instantiating the executor in a with block means we don't have to bother with shutting it down manually.

    # currently set to one as having it multithreaded seems to intermittently cause MTBLC folders to be given
    # incorrect id's
    max_workers = 2
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executioner:
        executioner.map(arg_wrapper, compoundList)

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-l', '--launch_directory', action=readable_dir, default = "" )
    parser.add_argument('-w', '--destination', action=readable_dir, help="Output directory", default=ebi_default_dest)
    parser.add_argument('-c', '--compound', help="- MetaboLights Compound Identifier", default="all")
    parser.add_argument('-s', '--stars', help="- Compounds import level", default="0")
    parser.add_argument('-f', '--ftp', action=readable_dir, default=ebi_ftp, help="FTP directory")
    parser.add_argument('-p', '--process', action=readable_dir, default="false", help="Use parallel threads")
    parser.add_argument('-d', '--debug', action='store_true', default=False, help='Setting debug to true will process only 20 compounds')
    args = parser.parse_args(arguments)

    # these variables being global stinks
    global workingDirectory
    global destinationDirectory
    global requestedCompound
    global ftp
    global importLevel

    debug = True # leave True to prevent concurrency (which causes mismatches currently)
    debug_incorrect_mtblc_dirs = False

    workingDirectory = args.launch_directory
    destinationDirectory = args.destination
    requestedCompound = args.compound.replace('"','')
    ftp = args.ftp
    importLevel = args.stars
    cron_debug = args.debug
    print(cron_debug)

    if(workingDirectory == ""):
        workingDirectory = os.getcwd()

    # log file configuration
    st = utils.getDateAndTime()
    randomInt = str(randint(1, 1000))
    logDirectory = workingDirectory + "/logs/" + st 
    if not os.path.exists(logDirectory):
        os.makedirs(logDirectory)
    logging.basicConfig(filename= logDirectory + "/log_" +randomInt +".log",level=logging.DEBUG)
    utils.init(logging)
    logging.info("-----------------------------------------------")
    logging.info('# Run started -' + utils.getDateAndTime())

    logging.info('Reading MetaboLights Study - Compound Mapping file')

    global mlSCMappingFile
    mlSCMappingFile = ftp + "mapping.json"

    logging.info('Reading Reactome data')

    global reactomeJSONFile
    reactomeJSONFile = ftp + "reactome.json"

    with open(reactomeJSONFile) as reactome_file:
        global reactomeData
        reactomeData = json.load(reactome_file)
    
    with open(mlSCMappingFile) as mapping_file:  
        global mlMapping
        mlMapping = json.load(mapping_file)

    if (requestedCompound != "all") :
        session = requests.Session()
        utils.fetchCompound(requestedCompound.strip(), workingDirectory, destinationDirectory, reactomeData, mlMapping, session)
    else:
        start_time = time.time()
        if debug_incorrect_mtblc_dirs:
            with open('bad_ids.txt') as bad_ids:
                list = ['MTBLC' + line.replace("\n","") for line in bad_ids.readlines()]
        else:
            list = utils.fetchMetaboLightsCompoundsList()
        if len(list) == 0:
            raise TimeoutError
        if debug is True:
            count = 0
            for compound in list:
                session = requests.Session()
                utils.fetchCompound(
                    compound.strip(), workingDirectory, destinationDirectory, reactomeData, mlMapping, session, file_lock)
                count += 1
                if cron_debug is True and count == 50000:
                    break
        else:
            get_all_compounds(list)
        logging.info('completed in {}'.format(str(time.time() - start_time)))

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
