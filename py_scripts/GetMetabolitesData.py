#!/usr/bin/env python

# Author: Venkata chandrasekhar Nainala 
# Version: 0.1.0
# Email: mailcs76@gmail.com / venkata@ebi.ac.uk

""" A python crawler to iterate MetaboLights Studies/Compounds and execute specified scripts
"""

import os
import sys
import json
import threading

import requests

import utils
import argparse
from os.path import basename
import subprocess

root = ""
destination = ""
request = ""
ftp = ""
metabolites = []
batch = 1
MetaboLightsUrl = "http://www.ebi.ac.uk/metabolights/webservice/"
MetaboLightsCompoundsList = MetaboLightsUrl + "compounds/list"

cal_default_dest = '/Users/cmartin/projects/fake_compound_dir/'
cal_ftp = '/Users/cmartin/projects/fakeftppub/'

ebi_default_dest = '/nfs/www-prod/web_hx2/cm/metabolights/prod/reference/'
ebi_ftp = '/net/isilonP/public/rw/homes/tc_cm01/from_ftp_pub/compounds/'

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

def fetchMetaboLightsCompoundsList():
    global metabolites, MetaboLightsUrl, MetaboLightsCompoundsList
    metabolites = requests.get(MetaboLightsCompoundsList).json()['content']
    #metabolites = json.loads(response.read())['content']

thread_storage = threading.local()

def get_session():
    if not hasattr(thread_storage, 'session'):
        thread_storage.session = requests.Session()

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-l', '--launch_directory', action=readable_dir, default="")
    parser.add_argument('-w', '--destination', action=readable_dir, default=ebi_default_dest, help="Output directory")
    parser.add_argument('-f', '--ftp', action=readable_dir, default=ebi_ftp, help="FTP directory")
    parser.add_argument('-d', '--debug', action="store_true", default=False, help="Setting debug to true will process only 20 studies")
    args = parser.parse_args(arguments)

    global metabolites, root, destination, request, batch

    # Reading lauching directory and log file details
    ftp = args.ftp
    root = args.launch_directory
    destination = args.destination
    debug = args.debug

    if(root == ""):
        root = os.getcwd()

    mlSCMappingFile     = ftp + "mapping.json"
    reactomeJSONFile    = ftp + "reactome.json"

    utils.generateMLStudyCompoundMappingFile(mlSCMappingFile, debug)
    utils.getReactomeData(reactomeJSONFile)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
