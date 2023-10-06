import concurrent.futures
import errno
import json
import logging
import math
import os
import xml.etree.ElementTree as ET

from pydantic import BaseModel
from requests import Session
from requests.exceptions import SSLError, ConnectionError, HTTPError, Timeout
from functools import wraps


class RuntimeFlags(BaseModel):
    spectra: bool = True
    citations: bool = True
    wikipathways: bool = True
    kegg: bool = True
    cactus: bool = True
    rhea: bool = True
    verbose_logging: bool = False
    timeout: int = 900
    mapping: dict = {
        'ms_from_mona_wrapper': 'spectra',
        'citation_wrapper': 'citations',
        'wikipathways_wrapper': 'wikipathways',
        'kegg_wrapper': 'kegg',
        'cactus_wrapper': 'cactus',
        'reactions_wrapper': 'rhea'
    }


class CompoundBuilderObjs(BaseModel):
    """
    A collection of various keys and mappings that the builder refers to. Previously these were scattered global
    variables, or defined in-method. Hopefully collecting them all here will make understanding this script easier than
     understanding v1. If these various objects don't make sense in isolation, follow the script flow and they might
     make some more sense. Also, view a MTBLS1234.json file and that also might shed some light.
    """
    chebi_ns_map: dict = {
        "envelop": "http://schemas.xmlsoap.org/soap/envelope/",
        "chebi": "{http://www.ebi.ac.uk/webservices/chebi}"
    }
    chebi_basic_keys: list = [
        'definition', 'smiles', 'inchi', 'inchiKey', 'charge', 'mass', 'monoisotopicMass', 'chebiAsciiName',
    ]
    chebi_adv_keys_map: dict = {
        'Synonyms': 'synonyms',
        'IupacNames': 'iupac_names',
        'Formulae': 'formulae',
        'Citations': 'citations',
        'DatabaseLinks': 'database_links',
        'CompoundOrigins': 'compound_origins',
        'Species': 'species'
    }

    chebi_citation_keys = ['source', 'type', 'value']
    chebi_citation_keys_map: dict = {'source': 'source', 'type': 'type', 'value': 'data'}
    epmc_citation_keys_map: dict = {
        'title': 'title',
        'doi': 'doi',
        'abstract': 'abstractText',
        'author': 'authorString',
    }

    chebi_species_keys: list = ['SpeciesAccession', 'SourceType', 'SourceAccession']
    chebi_species_via_mapping_file_map: dict = {
        'Species': None,
        'SpeciesAccession': 'study',
        'MAFEntry': 'mafEntry',
        'Assay': 'assay'
    }

    # recognise the below is confusing, it is a quirk from the old script, where the key of the database link dict
    # doesn't match the xml.find search term. So we have a mapping of `database link dict key: search term`.
    chebi_database_link_map: dict = {
        'source': 'type',
        'value': 'data'
    }

    # map our metabolights compound dict keys to chebi compound dict keys
    # this might also seem a little confusing, as we have two dicts with similar keys. I have lifted a lot of this
    # from version 1, and hope to remove it in future updates.
    ml_compound_chebi_compound_map: dict = {
        "name": "chebiAsciiName",
        "definition": "definition",
        "iupacNames": "IupacNames",
        "smiles": "smiles",
        "inchi": "inchi",
        "inchiKey": "inchiKey",
        "charge": "charge",
        "averagemass": "mass",
        "exactmass": "monoisotopicMass",
        "formula": "Formulae",
        "species": "Species",
        "synonyms": "Synonyms"
    }

    # specify what kind of empty value to give for a given key if the type associated with that key is not a string
    ml_compound_absent_value_type_map: dict = {
        "iupacNames": [],
        "species": [],
        "synonyms": []
    }

    reactions_keys: dict = {
        "name": "equation",
        "id": "id",
        "biopax2": "biopax2",
        "cmlreact": "cmlreact"
    }


class MtblsWsUrls(BaseModel):
    metabolights_ws_url: str = "http://www.ebi.ac.uk/metabolights/webservice/"
    metabolights_ws_study_url: str = f'{metabolights_ws_url}study/'
    metabolights_ws_studies_list: str = f'{metabolights_ws_study_url}study/list'
    metabolights_ws_compounds_url: str = f'{metabolights_ws_url}compounds/'
    metabolights_ws_compounds_list: str = f'{metabolights_ws_compounds_url}list'


class KeggUrls(BaseModel):
    kegg_api: str = "http://rest.kegg.jp/conv/compound/chebi:"
    kegg_pathways_list_api: str = "http://rest.kegg.jp/link/pathway/"
    kegg_pathway_api: str = "http://rest.kegg.jp/get/"


class MiscUrls(BaseModel):
    chebi_api: str = "https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId="
    cts_api: str = "http://cts.fiehnlab.ucdavis.edu/service/compound/"
    cactus_api = "https://cactus.nci.nih.gov/chemical/structure/"
    epmc_api: str = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="
    reactome_url: str = "http://www.reactome.org/download/current/ChEBI2Reactome.txt"
    rhea_api: str = "https://www.rhea-db.org/rhea/"
    wikipathways_api: str = "https://webservice.wikipathways.org/findPathwaysByXref?ids="
    mona_api: str = "http://mona.fiehnlab.ucdavis.edu/rest/spectra/search?query=compound.metaData=q=%27name==\%22InChIKey\%22%20and%20value==\%22"
    new_mona_api: str = "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search?endpoint=search&query=exists(compound.metaData.name%3A%27InChIKey%27%20and%20compound.metaData.value%3A%27{0}%27)%20and%20exists((tags.text%3A%27LC-MS%27))&size=10"


class CompoundBuilderUrls(BaseModel):
    mtbls: MtblsWsUrls = MtblsWsUrls()
    kegg: KeggUrls = KeggUrls()
    misc_urls: MiscUrls = MiscUrls()


class CompoundBuilderConfig:
    objs: CompoundBuilderObjs = CompoundBuilderObjs()
    urls: CompoundBuilderUrls = CompoundBuilderUrls()
    rt_flags: RuntimeFlags = RuntimeFlags()


# wrote a decorator to save try/except-ing the same exception in every method
def xml_exception_angel(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ET.ParseError as e:
            print(f'XML parsing error occurred: {str(e)}')
        except AttributeError as e:
            print(f'Attribute error while calling .find on xml document: {str(func)} {str(e)}')

    return wrapper


# wrote a decorator to catch common http request exceptions
def http_exception_angel(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except KeyError as e:
            print(f'KeyError when processing http request : {str(e)}')
        except Timeout as e:
            print(f'Request timed out when processing http request : {str(e)}')
        except HTTPError as e:
            print(f'HTTP Error code received when processing http request : {str(e)}')
            print(f'args: {str(args)}')
        except SSLError as e:
            print(f'Secure Sockets Layer Error when processing http request. Is the endpoint HTTPS enabled?')
        except ConnectionError as e:
            logging.exception(str(e))
            print(f'args: {str(args)}')
            print(f'ConnectionError when processing http request, is server reachable?: {str(e)}')
        except json.decoder.JSONDecodeError as e:
            print(f'JSONDecode error in {str(func)} args:{args} kwargs{kwargs}: {str(e)}')

    return wrapper


# I think once the above configs are finalised i will split them out. Script is ok readable as is but want to spare my
# successors the pain.
########################################################################################################################
#                                                                                                                      #
#                            Multi-Stage Compound Reference Layer Building Script                                      #
#                                                                                                                      #
########################################################################################################################

"""
This script performs multiple stages of building up an MTBLC compound directory. When run for all MTBLC ID's, the 
overall output constitutes MetaboLights' compound reference layer. Written by callum martin in March 2023. If I am no 
longer at the EBI and you have some questions, feel free to contact me at callumceltic@gmail.com and I will do my best
to help.

Usage:
    python3 StartCompoundBuilder.py --ftp <path_to_reference_files> --destination <path_to_output_dir>
Cron Usage:
    bsub -u metabolights-dev -J "newcompounddirbuild" -q standard -R "rusage[mem=128000]" -M 128000 
    /nfs/production/odonovan/scripts/new_compound_bot.sh

Parameters:
    --ftp: str
        The path to the reference files: mapping.json and reactome.json. Mapping.json is a map of studies to compounds
        to species, so which studies are associated with which compounds, and which species are associated with those 
        compounds. Reactome.json is a cache of reactome pathway data for each chebi compound. Read more about reactome 
        data at https://reactome.org/ .
    --destination: str
        The path to the output directory where the processed compound directories will be saved. Mass Spec .json files 
        will also be saved to the same directory. 

Stages:
    Stage 1: Initialisation 
        This stage prepares several objects needed for the builder to run. These are:
            - CompoundBuilderConfig object, which contains all internal & external API endpoints, as well as keys and 
            maps used to translate API responses to our format.
            - Session object for quicker http requests, shared among threads.
            - Mapping dict, loaded into memory from the mapping.json file
            - Reactome dict, loaded into memory from the reactome.json file.
            - List of all MTBLC id's, retrieved from our legacy java webservice
            - Initialises the compound dict that the script builds up.
            
    Stage 2: Get Chebi Data
        This stage hits the Chebi API for this MTBLC compound ID. Our IDs (excluding the MTBLC portion) are derived 
        directly from chebi. The response from the chebi API is in xml, so we use xml.ET.find to extract out the 
        information we want. In get_chebi_data, the returned chebi dict is populated in two sub-stages. One populates
        'basic' information, where we can just copy substrings out of the xml response into the chebi dict. The second
        stage uses a class called ChebiPopulator, which handles more 'advanced' information extraction. In some cases 
        we need to iterate over a xml.find result, or we need to build up dictionary objects to represent more 
        complicated data models (such as citations or database links). 
        I have tried to use dict comprehensions where possible, with the reference keys stored in the 
        CompoundBuilderConfig. I find this to be readable and clean, compared to v1 of this script where we had endless 
        try / except blocks, one for each attempt to assign a single field to a dictionary.
        Once the advanced dict is ready, the basic and the advanced dict are returned in a new dict together using
         pythons dict unpacking operator like so {**dict1, **dict2}
         
    Stage 3: Merge chebi dict with our dict
        This stage merges the result of the previous stage with the compound dict we initialised in stage 1. This takes 
        place entirely within (an admittedly complex) dict comprehension. The reason for the chebi dict and our compound
        dict having different keys is that the method was rewritten fairly close to the original implementation, and I 
        didn't want to tamper with it too much. A future iteration of this script could have the chebi dict be populated
        in such a manner that it could just be merged directly with our compound dict, without the need for a big dict 
        comprehension and storing a bunch of keys and maps in config.
    
    Stage 4: Multi threaded external API calls
        This is the main change from version 1 of the script. Previously, calls to external API endpoints were made in 
        sequence. With this version of the script, each external API that we hit has a dedicated thread for doing so. 
        This is managed in part by pythons relatively new ThreadPoolExecutor class (read more about the 
        ThreadPoolExecutor class here: https://superfastpython.com/threadpoolexecutor-in-python/ ). The RuntimeFlags 
        config object dictates which external API threads will be enabled. If a flag is set to false, we insert an 
        'empty' dict into the results.
        
        The external API's that we currently hit are:
            - KEGG: https://www.kegg.jp/ for pathway information.
            - Wikipathways: https://www.wikipathways.org/ for pathway information.
            - EuropePMC: https://europepmc.org/ for citation information.
            - Cactus: https://cactus.nci.nih.gov/chemical/structure for compound structure information.
            - Rhea: https://www.rhea-db.org/ for reactions information.
            - MoNa: https://mona.fiehnlab.ucdavis.edu/ for spectra data.
        
        To add to this list, you will need to implement an API function in ExternalAPIHitter, a wrapper for that API 
        function (more detail on wrappers below), and configure the inputs for that function in the ataronchronon,
         function and add the new wrapper function the list functions to be given to a thread, also in ataronchronon. 
        The actual API endpoint should be held within a config object, like the others, so we have all endpoints in one 
        place for reference. 
        
        In the ataronchronon function, the inputs for each API thread are collected into a tuple ( as ThreadPoolExecutor
        threads to my understanding only take a single input), and those tuples are collected into a list. We also 
        collect each 'wrapper' function into a list , where each 'wrapper' function simply unpacks the tuple of inputs, 
        and calls the actual API method with those inputs. Each wrapper/actual function combination can be found in the 
        class ExternalAPIHitter. 
        Using these two lists, we create a list of 'futures' objects where each future is the eventual output of a given
         thread. Each thread returns a standard dict that follows the format {'name': 'spectra', 'results': {...}}. 
         Once each thread has completed, the list of these result dicts is returned.
    
    Stage 5: Sorting results from multi threaded process
        The results from the previous stage are somewhat raw, and need further processing and sorting before being added
        to the compound dict. To this end we have a class called ExternalAPIResultSorter. It takes in the results of the
        multithreaded process, and the in progress compound dict and iterates over each result dict. Using the name of 
        the result dict, it retrieves and calls the associated 'handling' function (one for each external API). The 
        handling function processes the results, and inserts them into the compound dict. Any empty results are logged 
        out and the relevant field in the dict set to an empty signifier. 'Flags' are also set here, indicating the 
        presence of a type of data. These 'flags' are used in the UI to decide whether to try and render a certain 
        portion of the compound data.
    
    Stage 6: Final odds and ends, and saving the dict to file
        There are some bits of data outstanding from the previous version that I couldn't integrate elsewhere. They are:
            - Processing reactome pathways from the reactome.json file. This is done in a single function called 
            get_reactome_data, and that function just tries to retrieve the slice corresponding to the current compound 
            ID.
            - Processing NMR data for this compound, using the results of our Compounds java webservice call. The legacy
            java web project also holds some compound information (aside from the list of ID's), which we use here 
            solely for populating the NMR field of the compound dict. We do this in a function called get_nmr
        Once the above functions have run, we set the flags that were not accounted for in the ExternalAPIResultSorter
        .sort function. These flags are 'hasPathways', 'hasNMR' and 'hasSpecies'. Again, these flags drive the UI by 
        telling it whether or not to attempt to load a particular portion of compound data.
        At this point, the compound is saved as a .json file in the destination directory in its own unique subdirectory
        . The dict itself is also returned, but nothing is done with it in the script that calls it 
        ( StartCompoundBuilder.py ), but is useful if you want to do some debugging, so I have left it as is.

Further Notes:
    Worth checking out the different ExternalAPIHitter methods to see what they do. One, get_ms_from_mona, includes file
    I/O operations. As a result this one method is a significant drag on the multithreaded process. Future versions of 
    this script could do with splitting out this portion entirely.
    
    There are some custom decorators I have written to save writing the same try except blocks over and over. While it
    obscures the root cause a little, you can see for which compound and what method the exception was thrown on. Since 
    we are processing so much data, I decided this level of exception handling was ok.
    
    Some tests would be good, could even have one that built 50 or so real compound directories locally for 
    scrutinisation. Not completed at the time due to other high priority tasks.
    
    The venv for this script must be maintained manually. Any new dependencies will require you to login to codon and
    install the requirements.txt again (obviously make your your new dependency is added to requirements.txt, and 
    transfer that to codon, to allow us to keep track at least somewhat).
    
    If you want to make this script single threaded for whatever reason, you can change the following line:
    `concurrent.futures.ThreadPoolExecutor(max_workers=x)` in ataronchron, and replace whatever the current value of x 
    is to 1.

Output:
    The script saves the compound and any intermediate data files to the specified output directory.

"""

def build(metabolights_id, ml_mapping, reactome_data, data_directory):
    """
    Entrypoint method for the script, to build an MTBLC compound directory.

    :param metabolights_id: the MTBLC12345 ID retrieved from the web service.
    :param ml_mapping: The mapping file, which associates studies to compound ids referenced in that study.
    :param reactome_data: Reactome data json file.
    :param data_directory: Directory to save the built compound subdirectory to.

    :return: N/A but saves the directory to the data directory.
    """
    config = CompoundBuilderConfig()
    session = Session()
    chebi_id = metabolights_id.replace("MTBLC", "").strip()

    # call our java webservice
    mtblcs = None
    try:
        mtblcs = session.get(f'{config.urls.mtbls.metabolights_ws_compounds_url}{metabolights_id}').json()['content']
    except json.JSONDecodeError as e:
        print(f'Error getting info from MTBLS webservice for compound {chebi_id}')
    if mtblcs is None:
        print(f'Exiting compound building process for compound {chebi_id}')
        return {}
    # init our compound dict, build the chebi dict
    compound_dict = _InternalUtils.initialize_compound_dict()
    chebi_dict = get_chebi_data(chebi_id, ml_mapping, config, session)
    if not chebi_dict:
        return compound_dict
    chebi_dict['id'] = chebi_id
    # do some updating
    compound_dict['id'] = metabolights_id

    # perform the gnarliest dict comprehension anyone has ever seen
    # essentially just copies values from the chebi compound dict to the mtbl dict, referring to
    # `ml_compound_chebi_compound_map` so it knows which key on the chebi dict matches which key on the mtbl dict.
    # Also, if a value is not present, and the type of that value is not a string, it refers to the
    # `ml_compound_absent_type_value` map, which specifies what kind of empty value to concatenate
    compound_dict.update(
        {
            key: chebi_dict.get(
                value, 'NA'
                if key not in config.objs.ml_compound_absent_value_type_map
                else config.objs.ml_compound_absent_value_type_map[key]
            )
            if chebi_dict.get(value) is not None
            else (print(f"{value} not assigned"), 'NA'
            if key not in config.objs.ml_compound_absent_value_type_map
            else config.objs.ml_compound_absent_value_type_map[key])[1]
            for key, value in config.objs.ml_compound_chebi_compound_map.items()
        }
    )

    # initialise the pathways dicts and the spectra lists
    compound_dict.update({
        'pathways': {'WikiPathways': {}, 'KEGGPathways': {}, 'ReactomePathways': {}},
        'spectra': {'NMR': [], 'MS': []}
    })

    mementos = ataronchronon(
        chebi_compound_dict=chebi_dict,
        config=config,
        session=session,
        mtbls_id=metabolights_id,
        data_directory=data_directory
    )

    apocrypha = ExternalAPIResultSorter(mementos)
    compound_dict = apocrypha.sort(compound_dict)

    # last bits of data that couldn't be integrated anywhere else added here
    compound_dict['pathways']['ReactomePathways'] = get_reactome_data(metabolights_id, reactome_data)
    compound_dict['spectra']['NMR'] = get_nmr(mtblcs['mc']['metSpectras']) if mtblcs else []

    # update NMR, species, pathways flags.
    if (len(compound_dict['pathways']['ReactomePathways']) + len(compound_dict['pathways']['KEGGPathways']) +
        len(compound_dict['pathways']['WikiPathways'])) > 0:
        compound_dict['flags']['hasPathways'] = 'true'
    if len(compound_dict['spectra']['NMR']) > 0:
        compound_dict['flags']['hasNMR'] = 'true'
    if compound_dict['species']:
        compound_dict['flags']['hasSpecies'] = 'true'

    if config.rt_flags.verbose_logging:
        print(f'___________________________ataronchoron results for {metabolights_id}_____________')
        for d in mementos:
            print(d.values())

    _FileHandler.save_json_file(f'{data_directory}/{metabolights_id}/{metabolights_id}_data.json', compound_dict)
    return compound_dict


def ataronchronon(chebi_compound_dict: dict, config: CompoundBuilderConfig, session: Session, mtbls_id: str,
                  data_directory: str):
    """
    Configure a ThreadPoolExecutor, and instruct each thread to execute a different external API related task.
    The threads for each task will only start if the corresponding flag is enabled. The RuntimeFlags object starts
    with all flags enabled by default.
    :param ml_compound_dict: The in-progress MTBLC compound dict.
    :param chebi_compound_dict: compound dict built from results of chebi API response earlier in compound building.
    :param config: CompoundBuilderConfig object.
    :param session: Session object, shared among threads, to make http calls.
    :param mtbls_id: the MTBLC12335 ID of the current compound.
    :param data_directory: Where the compound directory and spectra files will be saved.
    :return: list of result dicts.
    """
    # collect the inputs for each method wrapper wrapper into a tuple. We do this as threadpool executors can only take
    # one argument.
    citation_input = (chebi_compound_dict['Citations'], config, session)
    cactus_input = (config.urls.misc_urls.cactus_api, chebi_compound_dict['inchiKey'], session)
    reactions_input = (chebi_compound_dict, config.urls.misc_urls.rhea_api, config.objs, session)
    ms_from_mona_input = (mtbls_id, data_directory, chebi_compound_dict['inchiKey'], config, session)
    wiki_pathways_input = (chebi_compound_dict['inchiKey'], config, session)
    kegg_pathways_input = (chebi_compound_dict, config, session)

    input_list = [
        citation_input, cactus_input, reactions_input, ms_from_mona_input, wiki_pathways_input, kegg_pathways_input
    ]
    method_list = [ExternalAPIHitter.citation_wrapper, ExternalAPIHitter.cactus_wrapper,
                   ExternalAPIHitter.reactions_wrapper, ExternalAPIHitter.ms_from_mona_wrapper,
                   ExternalAPIHitter.wikipathways_wrapper, ExternalAPIHitter.kegg_wrapper]

    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ur_executor:
        # create a list of empty result dicts to make processing easier later
        the_duds = [
            {'name': config.rt_flags.mapping[_InternalUtils.extract_name_from_function(method)], 'results': None}
            for method in method_list
            if _InternalUtils.flag_is_enabled(config.rt_flags, method) is False
        ]

        # create a list of futures object, where each future is a thread corresponding to an enabled runtime flag.
        the_futures = [
            ur_executor.submit(method, args) for method, args in zip(method_list, input_list)
            if _InternalUtils.flag_is_enabled(config.rt_flags, method)
        ]

        # wait for and collect the results of each individual thread
        the_results = [future.result() for future
                       in concurrent.futures.as_completed(the_futures, config.rt_flags.timeout)]

        # add the duds to the results
        the_results.extend(the_duds)
        return the_results


@xml_exception_angel
def get_chebi_data(id, ml_mapping, config, session: Session) -> dict:
    """
    Hit the ChEBI API and parse the response using xml ElementTree. It then populates a 'basic' dict from the initial
    response (basic as all it does is make a single xml.find query per key, and the keys can be found in
    `CompoundBuilderObjs.chebi_basic_keys`). A ChebiPopulator instance is then used to build an 'advanced' chebi dict
    (advanced as we need a bespoke method to extract the necessary information for each  key, and the keys can be found
    in `CompoundBuilderObjs.chebi_adv_keys_map`). The two dicts are combined in the return statement.

    :param id: chebi ID that we want all their information on.
    :param ml_mapping: mapping object that links studies to chebi ID's and species.
    :param config: CompoundBuilderConfig instance.
    :param ml_mapping: The mapping file, which associates studies to compound ids referenced in that study.

    :return: dict representing the information held in chebi on a particular compound.
    """
    """This request needs to be made first, as the rest of the other API calls depend on it's response."""

    chebi_response = session.get(f'{config.urls.misc_urls.chebi_api}{id}').content
    root = ET.fromstring(chebi_response) \
        .find("envelop:Body", namespaces=config.objs.chebi_ns_map) \
        .find("{https://www.ebi.ac.uk/webservices/chebi}getCompleteEntityResponse") \
        .find("{https://www.ebi.ac.uk/webservices/chebi}return")

    # log out the ID so we know it's a valid 'un
    _InternalUtils.gimme_line()
    print(root.find("{https://www.ebi.ac.uk/webservices/chebi}chebiId"))

    # generate the 'basic' dict.
    chebi_basic_dict = {
        key: _InternalUtils.get_val(root, key)
        for key in config.objs.chebi_basic_keys
    }
    chebi_basic_dict['id'] = id

    # init the ChebiPopulator class, and chain call its methods
    chebi_advanced_populator = ChebiPopulator(root, config)
    chebi_advanced_populator \
        .get_synonyms() \
        .get_iupac_names() \
        .get_formulae() \
        .get_citations() \
        .get_database_links() \
        .get_species_via_compound_origins() \
        .get_species_via_compound_mapping(ml_mapping, chebi_basic_dict['id'])

    chebi_advanced_dict = {
        key: chebi_advanced_populator.__getattribute__(value)
        for key, value in config.objs.chebi_adv_keys_map.items()
    }

    # mash the two dicts together
    return {**chebi_basic_dict, **chebi_advanced_dict}


def get_reactome_data(mtblc_compound_id: str, reactome_data: dict) -> dict:
    """
    Go over the pathways in the reactome file that correspond to this particular compound. Return the results
    as a dict.
    :param mtblc_compound_id: Compound ID to slice the reactome data with.
    :param reactome_data: Reactome data as a dict.
    :return: Reactome pathways for this compound, as a dict.
    """
    temp_reactome_pathways = reactome_data[mtblc_compound_id] if mtblc_compound_id in reactome_data else []
    reactome_pathways = {}
    try:
        for pathway in temp_reactome_pathways:
            temp_pathway = {'name': pathway['pathway'], 'pathwayId': pathway['pathwayId'],
                            'url': pathway['reactomeUrl'],
                            'reactomeId': pathway['reactomeId']}
            if pathway['species'] not in reactome_pathways:
                reactome_pathways[pathway['species']] = [temp_pathway]
            else:
                reactome_pathways[pathway['species']].append(temp_pathway)
    except KeyError as e:
        print(f'Error populating dict for {mtblc_compound_id}: {str(e)}')
    finally:
        return reactome_pathways


def get_nmr(spectra) -> list:
    """
    Iterate over the java webservice result for this compound, and populate a dict for each spectra in the results.
    Any KeyErrors are logged and the offending spectra skipped.
    :param spectra: Spectra info from java webservice.
    :return: List of NMR spectra dicts.
    """
    nmr = []
    for spec in spectra:
        try:
            if spec["spectraType"] == "NMR":
                temp_spec = {'name': spec['name'], 'id': str(spec['id']),
                             'url': f"http://www.ebi.ac.uk/metabolights/webservice/compounds/spectra/{str(spec['id'])}/json",
                             'path': spec['pathToJsonSpectra'], 'type': spec['spectraType'],
                             'attributes': [
                                 {'attributeName': attr['attributeDefinition']['name'],
                                  'attributeDescription': attr['attributeDefinition']['name'],
                                  'attributeValue': attr['value']} for attr in spec['attributes']
                             ]}
                nmr.append(temp_spec)
        except KeyError as e:
            print(f'Error populating dict for {spec["id"]}: {str(e)}')
            continue
    return nmr


class ChebiPopulator:

    def __init__(self, root, config: CompoundBuilderConfig):
        self.root = root
        self.config = config

        self.synonyms = []
        self.iupac_names = []
        self.formulae = None
        self.citations = []
        self.database_links = []
        self.compound_origins = []
        self.species = {}

    @xml_exception_angel
    def get_synonyms(self):
        """
        Search the chebi compound xml root for Synonyms, and then iterate over the results, further searching each
        result for a data tag, and then appending the contents of that data tag to the instances `synonyms` list.

        :return: Self to enable chain calling.
        """
        for synonym in self.root.findall("{https://www.ebi.ac.uk/webservices/chebi}Synonyms"):
            self.synonyms.append(synonym.find("{https://www.ebi.ac.uk/webservices/chebi}data").text)
        return self

    @xml_exception_angel
    def get_iupac_names(self):
        """
        Search the chebi compound xml root for Iupac Names, and then iterate over the results, further searching each
        result for a data tag, and then appending the contents of that data tag to the instances `iupac_names` list.

        :return: Self to enable chain calling.
        """
        for iupac_name in self.root.findall("{https://www.ebi.ac.uk/webservices/chebi}IupacNames"):
            self.iupac_names.append(iupac_name.find("{https://www.ebi.ac.uk/webservices/chebi}data").text)
        return self

    @xml_exception_angel
    def get_formulae(self):
        """
        Search the chebi compound xml root for Formulae, and then further search for a data tag within Formulae, and
        set the instances `formulae` attribute to that value.

        :return: Self to enable chain calling.
        """
        if self.root:
            self.formulae = self.root \
                .find("{https://www.ebi.ac.uk/webservices/chebi}Formulae") \
                .find("{https://www.ebi.ac.uk/webservices/chebi}data").text

        return self

    @xml_exception_angel
    def get_citations(self):
        """
        Search the chebi compound xml root for Citations, then iterate over the results, creating a new citation dict
        per result, then appending that dict to the instances `citations` list.
        The keys for the citation dict can be found in the CompoundBuilderObjs class definition.

        :return: Self to enable chain calling.
        """
        for citation in self.root.findall("{https://www.ebi.ac.uk/webservices/chebi}Citations"):
            citation_dict = {
                key: citation.find("{https://www.ebi.ac.uk/webservices/chebi}" + value).text
                if citation.find("{https://www.ebi.ac.uk/webservices/chebi}" + value) is not None else 'N/A'
                for key, value in self.config.objs.chebi_citation_keys_map.items()
            }
            self.citations.append(citation_dict)
        return self

    @xml_exception_angel
    def get_database_links(self):
        """
        Search the chebi compound xml root for Database Links, then iterate over the results, creating a new database
        link dict per result, then appending that dict to the instances `database_links` list.
        The keys for the database link dict can be found in the CompoundBuilderObjs class definition.

        :return: Self to enable chain calling.
        """
        for database_link in self.root.findall("{https://www.ebi.ac.uk/webservices/chebi}DatabaseLinks"):
            database_link_dict = {
                key: database_link.find(
                    "{https://www.ebi.ac.uk/webservices/chebi}" + self.config.objs.chebi_database_link_map[key]
                ).text
                for key in self.config.objs.chebi_citation_keys
                if key is not 'type'
            }
            self.database_links.append(database_link_dict)
        return self

    @xml_exception_angel
    def get_species_via_compound_origins(self):
        """
        Search the chebi compound xml root for Compound Origin species, then iterate over the results, creating a new
        species dict per result, then appending that dict to a list within the instance's `species` dict. That list will
        contain each result dict from this method and `get_species_via_compound_mapping`.
        The keys for the species dict can be found in the CompoundBuilderObjs class definition.

        :return: Self to enable chain calling.
        """
        for origin in self.root.findall("{https://www.ebi.ac.uk/webservices/chebi}CompoundOrigins"):
            chebi_species = origin.find("{https://www.ebi.ac.uk/webservices/chebi}speciesText").text.lower()
            if chebi_species not in self.species:
                self.species[chebi_species] = []
            origin_dict = {
                key: origin.find("{https://www.ebi.ac.uk/webservices/chebi}" +
                                 f"{_InternalUtils.pascal_case(key) if key == 'SpeciesAccession' else key}"
                                 ).text if origin.find("{https://www.ebi.ac.uk/webservices/chebi}" +
                                                       f"{_InternalUtils.pascal_case(key) if key == 'SpeciesAccession' else key}"
                                                       ) is not None else "N/A"
                for key in self.config.objs.chebi_species_keys
            }
            self.species[chebi_species].append(origin_dict)
        return self

    @xml_exception_angel
    def get_species_via_compound_mapping(self, mapping: dict, id: str):
        """
        I may rewrite this, passing the giant object to this method feels a bit wrong. This could be a crunch point also,
        as python will have to search the giant mapping object each time. Can't conceive a different way of getting the
        information out of that mapping file currently.
        Retrieves the species from the study mapping file, and populates a bunch of dicts, one per species, using info
        from the mapping file.

        :param mapping: Big study-compound-species mapping object.
        :param id: CHeBI compound ID
        :return: Self to enable chain calling.
        """
        if id in mapping:
            study_species = mapping['compoundMapping'][id]
            for study_s in study_species:
                temp_study_species = str(study_s['species']).lower()
                if temp_study_species not in self.species:
                    self.species[temp_study_species] = []
                origin_dict = {
                    key: study_s[self.config.objs.chebi_species_via_mapping_file_map[key]]
                    if key is not 'Species' else temp_study_species
                    for key in self.config.objs.chebi_species_via_mapping_file_map.keys()
                }
                self.species[temp_study_species].append(origin_dict)
        return self


class ExternalAPIHitter:
    """
    The methods in this class come in twos:
    - method_wrapper(tuple_of_inputs)
    - actual_method(and, its, inputs)

    The method wrapper is what gets passed to the ThreadPoolExecutor, and methods that get passed to the TPE can only
    take one argument, so we package all argument for a given method into a tuple, and the unpack it using the `*`
    operator when we call the method proper.
    """

    @staticmethod
    def citation_wrapper(citation_tuple) -> dict:
        return {'results': ExternalAPIHitter.get_citations(*citation_tuple), 'name': 'citations'}

    @staticmethod
    @http_exception_angel
    def get_citations(citations, config: CompoundBuilderConfig, session: Session) -> list:
        """
        For each citation for a given compound, hit the europePMC API, get json format of the result, and update the
        existing
        :param citations: citations that we have from the chebi response from earlier in the compound building process.
        :param config: CompoundBuilderConfig object, used to get dictionary keys and api endpoints.
        :param session: Session object to make http calls.
        :return: list of updated citations
        """
        print("Attempting to get data from europePMC API.")
        epmc_list = []
        for citation in citations:
            print(
                f'attempting to hit {config.urls.misc_urls.epmc_api}{str(citation["value"])}&format=json&resultType=core&cursorMark=*&pageSize=25')
            try:
                citation_epmc_data = session.get(
                    f'{config.urls.misc_urls.epmc_api}{str(citation["value"])}&format=json&resultType=core').json()[
                    'resultList']['result'][0]
            except json.decoder.JSONDecodeError as e:
                print(f'No response for individual citation {str(citation["value"])}:{str(e)}')
                continue
            except IndexError as e:
                print(f'No response for individual citation {str(citation["value"])}:{str(e)}')
                continue
            citation.update(
                {key: citation_epmc_data[value] if value in citation_epmc_data else 'NA'
                 for key, value in config.objs.epmc_citation_keys_map.items()}
            )
            epmc_list.append(citation)
        return epmc_list

    @staticmethod
    def cactus_wrapper(cactus_tuple) -> dict:
        return {'results': ExternalAPIHitter.get_cactus_structure(*cactus_tuple), 'name': 'cactus'}

    @staticmethod
    @http_exception_angel
    def get_cactus_structure(cactus_api: str, inchi_key: str, session: Session) -> str:
        """
        Hit the cactus API, and return the text of the result.
        :param cactus_api: cactus api endpoint to hit.
        :param inchi_key: inchi_key associated with the current compound.
        :param session: Sesion object used to make http call.
        :return: string of cactus API response.
        """
        print('Attempting to get data from cactus API.')
        return session.get(f'{cactus_api}{inchi_key}/sdf').text

    @staticmethod
    def reactions_wrapper(reactions_tuple) -> dict:
        return {'results': ExternalAPIHitter.get_reactions(*reactions_tuple), 'name': 'reactions'}

    @staticmethod
    @http_exception_angel
    def get_reactions(chebi_compound_dict, rhea_api, conf_objs: CompoundBuilderObjs, session: Session) -> list:
        """
        Ping the RHEA API with the chebi ID for the current compound as a query parameter. Then, for each result, parse
        that result into a dict, and include that dict in the list of reactions. Note that the dict comprehension that
        builds a single dict, and the list comprehension which builds the list of reaction dicts take place at the same
        time.
        In config.reaction_keys.items(), there are two keys, biopax2 and cmlreact that ultimately go unused. This is
        because version 1 of the script hit a long outdated version of the rhea api, that would also given chemical
        markup language (cml) output, and biopax2 information in it's response. The current rhea api does not offer this
        information, but I have left the keys in for posterity, in case we discover another way to get the missing
        information.
        :param chebi_compound_dict: compound dict built from results of chebi API response earlier in compound building.
        :param rhea_api: Endpoint for the swiss bioinformatics institute's rhea api
        :param conf_objs: CompoundBuilderObjs object used to get reaction dict keys.
        :param session: Session object to make http call.
        :return: list of reaction dicts.
        """
        print('Attempting to get data from rhea API.')
        query = "?query="
        columns = "&columns=rhea-id,equation,chebi-id"
        format = "&format=json"
        limit = "&limit=10"
        rhea_data = session.get(f'{rhea_api}{query}{chebi_compound_dict["id"]}{columns}{format}{limit}').json()
        print(f'rhea data for chebi id {chebi_compound_dict["id"]} : {rhea_data["results"]}')
        reactions = [
            {
                key: result[value] if value in result else ""
                for key, value in conf_objs.reactions_keys.items()
            } for result in rhea_data['results']
        ]
        return reactions

    @staticmethod
    def ms_from_mona_wrapper(spectra_tuple) -> dict:
        return {'results': ExternalAPIHitter.get_ms_from_mona(*spectra_tuple), 'name': 'spectra'}

    @staticmethod
    @http_exception_angel
    def get_ms_from_mona(
            mtbls_id: str, dest: str, inchi_key: str, config: CompoundBuilderConfig, session: Session) -> list:
        """
        Ping the MoNA API with the inchikey for a given compound as a query parameter. Then, for each result, parse that
        result into a spectra dict, and pass that dict to `_FileHandler.save_spectra` to further process the spectral
        data and save it as a .json file.
        :param mtbls_id: The MTBLC accession number associated with the MTBLC compound we are building.
        :param dest: The parent compound reference directory.
        :param inchi_key: inchi_key associated with a given compound.
        :param config: CompoundBuilderConfig object, mona_api endpoint extracted from within.
        :param session: Session object to make http call.
        :return: list of spectra objects representing a spectrum.
        """
        print(
            f"Attempting to get spectral data from MoNa at {f'{config.urls.misc_urls.new_mona_api.format(inchi_key)}'}")
        ml_spectrum = []
        mona_api_url  = f'{config.urls.misc_urls.new_mona_api.format(inchi_key)}'
        result = session.get(mona_api_url).json()
        for spectra in result:
            ml_spectra = {'splash': spectra['splash'], 'type': 'MS', 'name': str(spectra['id']),
                          'url': f'/metabolights/webservice/beta/spectra/{mtbls_id}/{str(spectra["id"])}'}
            temp_submitter = spectra['submitter']
            ml_spectra['submitter'] = f"{str(temp_submitter['firstName'])}  {str(temp_submitter['lastName'])} ; " \
                                      f"{str(temp_submitter['emailAddress'])} ; {str(temp_submitter['institution'])}"
            ml_spectra['attributes'] = []
            for metadata in spectra['metaData']:
                if not metadata['computed']:
                    temp_attribute = {'attributeName': metadata['name'], 'attributeValue': metadata['value'],
                                      'attributeDescription': ""}
                    ml_spectra['attributes'].append(temp_attribute)
            ml_spectrum.append(ml_spectra)
            _FileHandler.save_spectra(str(spectra['id']), spectra['spectrum'], mtbls_id, dest)
        return ml_spectrum

    @staticmethod
    def wikipathways_wrapper(pathway_tuple) -> dict:
        return {'results': ExternalAPIHitter.get_wikipathways(*pathway_tuple), 'name': 'wikipathways'}

    @staticmethod
    @http_exception_angel
    def get_wikipathways(inchi_key: str, config: CompoundBuilderConfig, session: Session) -> dict:
        """
        Hit the wikipathways API, and for each result / pathway, parse it into a new dict and append that pathway to a
        parent `final_pathways` object. The final object is the representation of all pathways for a given compound,
        sorted by species.
        :param inchi_key: inchi_key of a particular compound to use as a query param.
        :param config: CompoundBuilderConfig object to pull out the wikipathways url.
        :param session: Session object to make http call.
        :return: dict object representing pathways for particular compound.
        """
        format_params = "&codes=Ik&format=json"
        print(
            f'Attempting to retrieve wikipathways data from {config.urls.misc_urls.wikipathways_api}{inchi_key}{format_params}')
        final_pathways = {}
        wikipathways = session.get(
            f'{config.urls.misc_urls.wikipathways_api}{inchi_key}{format_params}').json()['result']

        for pathway in wikipathways:
            if pathway['species'] not in final_pathways:
                final_pathways[pathway['species']] = []
            pathway_dict = {'id': pathway['id'], 'url': pathway['url'], 'name': pathway['name']}
            if pathway_dict not in final_pathways[pathway['species']]:
                final_pathways[pathway['species']].append(pathway_dict)

        return final_pathways

    @staticmethod
    def kegg_wrapper(kegg_tuple) -> dict:
        return {'results': ExternalAPIHitter.get_kegg_pathways(*kegg_tuple), 'name': 'kegg_pathways'}

    @staticmethod
    @http_exception_angel
    def get_kegg_pathways(chebi_compound_dict: dict, config: CompoundBuilderConfig, session: Session) -> list:
        """
        Hit the kegg API using the compounds chebi ID as a query parameter. Then, hit kegg's pathway list API using the
        kegg id that we got from the previous request. Then, for each line in the pathways list response, hit kegg's
        individual pathway API for that line, and then parse the response into a pathway dict.
        :param chebi_compound_dict: compound dict built from results of chebi API response earlier in compound building.
        :param config: CompoundBuilderConfig object, used to pull out kegg api endpoints.
        :param session: Session object to make http calls.
        :return: list of kegg pathway objects.
        """
        print(f'Attempting to retrieve KEGG data from {config.urls.kegg.kegg_api}{chebi_compound_dict["id"].lower()}')
        final_kegg_pathways = []
        kegg_id_q_r = session.get(f'{config.urls.kegg.kegg_api}{chebi_compound_dict["id"].lower()}')
        kegg_id = None
        try:
            kegg_id = kegg_id_q_r.text.split("\t")[1].strip()
        except IndexError as e:
            print(f'Unable to get a corresponding CPD number for chebi compound {chebi_compound_dict["id"]}: {str(e)}')
        if kegg_id is None:
            return final_kegg_pathways

        print(f'Attempting step 2 of getting kegg data with url {config.urls.kegg.kegg_pathways_list_api}{kegg_id}')
        pathways_data = session.get(f'{config.urls.kegg.kegg_pathways_list_api}{kegg_id}').text
        for line in pathways_data.strip().split("\n"):
            if line == "":
                continue
            try:
                pathway_id = line.split("\t")[1].strip()
            except IndexError as e:
                print(f'Couldnt get pathway id due to index error when parsing pathways response: {str(e)}')
                continue
            pathway_data = session.get(f'{config.urls.kegg.kegg_pathway_api}{pathway_id}').text
            pathway_dict = {'id': pathway_id}
            for pline in pathway_data.strip().split("\n"):
                if "NAME" in pline:
                    pathway_dict['name'] = pline.replace("NAME", "").strip()
                elif "KO_PATHWAY" in pline:
                    pathway_dict['KO_PATHWAYS'] = pline.replace('KO_PATHWAYS', '').strip()
                elif "DESCRIPTION" in pline:
                    pathway_dict['description'] = pline.replace('DESCRIPTION', "").strip()
            final_kegg_pathways.append(pathway_dict)
        return final_kegg_pathways


class ExternalAPIResultSorter:

    def __init__(self, mementos):
        self.mementos = mementos

    def __getattribute__(self, name):
        """
        Instead of just returning the attribute as normal, we want to return a function object that we can call.
        :param name: Name of the function we want to call.
        :return: Lamdba function that calls the requested handling_function.
        """
        try:
            return super().__getattribute__(name)
        except AttributeError:
            handling_function = super().__getattribute__(name)
            return lambda x, y: handling_function(x, y)

    def sort(self, metabolights_dict: dict) -> dict:
        """
        Entry method for the ExternalApiResultsSorter. Iterates over each memento returned by the multithreaded
        `ataronchronon` process, and retrieves and calls the handling method for that memento. The handling method
        updates the metabolights_dict with the results of the memento (with some kind of formatted null value if no
        results are presented
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict updated with all results from the `ataronchronon` multithreaded process.
        """
        for memento in self.mementos:
            metabolights_dict = self.__getattribute__(f'handle_{memento["name"]}')(memento, metabolights_dict)
        return metabolights_dict

    def handle_cactus(self, cactus_memento, metabolights_dict: dict) -> dict:
        """
        Check if the cactus_memento has any empty result signifiers, and if it don't, assign the results
        to the metabolights_dict's structure field.
        :param cactus_memento: Output from the cactus thread of the ataronchronon process.
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict with the structure field updated.
        """
        if cactus_memento['results'] is None or cactus_memento['results'] == []:
            metabolights_dict['structure'] = 'NA'
            print(f'Compound Error {metabolights_dict["id"]} Structure not assigned.')
            return metabolights_dict

        metabolights_dict['structure'] = cactus_memento['results']
        return metabolights_dict

    def handle_citations(self, citations_memento, metabolights_dict: dict) -> dict:
        """
        Check if the citations_memento has any empty result signifiers, and if it don't, assign the results
        to the metabolights_dict's citations field. Also set the citations flag accordingly
        :param citations_memento: Output from the citations thread of the ataronchronon process.
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict with the citations field and flag updated.
        """
        if citations_memento['results'] is None or citations_memento['results'] == []:
            metabolights_dict['citations'] = []
            metabolights_dict['flags']['hasLiterature'] = 'false'
            return metabolights_dict

        metabolights_dict['citations'] = citations_memento['results']
        metabolights_dict['flags']['hasLiterature'] = 'true'
        return metabolights_dict

    def handle_spectra(self, spectra_memento, metabolights_dict: dict) -> dict:
        """
        Check if the spectra memento has any empty result signifiers, and if it don't, assign the results
        to the metabolights_dict's spectra['MS'] field. Also set the hasMS flag accordingly
        :param spectra_memento: Output from the spectra thread of the ataronchronon process.
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict with the spectra['MS'] field and flag updated.
        """
        if spectra_memento['results'] is None or spectra_memento['results'] == []:
            print(f'No MoNa info available for {metabolights_dict["id"]}')
            metabolights_dict['flags']['hasMS'] = 'false'
            return metabolights_dict

        metabolights_dict['spectra']['MS'] = spectra_memento['results']
        metabolights_dict['flags']['MS'] = 'true'
        return metabolights_dict

    def handle_kegg_pathways(self, kegg_memento, metabolights_dict: dict) -> dict:
        """
        Check if the kegg memento has any empty result signifiers, and if it don't, assign the results
        to the metabolight_dict's pathways['KEGGPathways'] field.
        :param kegg_memento: Output from the kegg thread of the ataronchronon process.
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict with the pathways['KEGGPathways'] field updated.
        """
        if kegg_memento['results'] is None or kegg_memento['results'] == {}:
            print(f'No KEGG info for {metabolights_dict["id"]}')
            return metabolights_dict
        metabolights_dict['pathways']['KEGGPathways'] = kegg_memento['results']
        return metabolights_dict

    def handle_wikipathways(self, wiki_memento, metabolights_dict: dict) -> dict:
        """
        Check if the wikipathways memento has any empty result signifiers, and if it don't, assign the results
        to the metabolights_dict's pathways['WikiPathways'] field.
        :param wiki_memento: Output from the wikipathways thread of the ataronchronon process.
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict with the pathways['WikiPathways'] field updated.
        """
        if wiki_memento['results'] is None or wiki_memento['results'] == {}:
            print(f'No WikiPathways info for {metabolights_dict["id"]}')
            return metabolights_dict
        metabolights_dict['pathways']['WikiPathways'] = wiki_memento['results']
        return metabolights_dict

    def handle_reactions(self, reactions_memento, metabolights_dict: dict) -> dict:
        """
        Check if the reactions memento has any empty result signifiers, and if it don't, assign the results
        to the metabolights_dict's reactions field. Also set the hasReactions flag accordingly.
        :param reactions_memento: Output from the reactions thread of the ataronchronon process.
        :param metabolights_dict: The in progress metabolights compound dict.
        :return: The metabolights_dict with the reactions field and flag updated.
        """
        if reactions_memento['results'] is None or reactions_memento == []:
            print(f'No Rhea info for {metabolights_dict["id"]}. Reactions not assigned.')
            metabolights_dict['flags']['hasReactions'] = 'false'
            return metabolights_dict

        metabolights_dict['reactions'] = reactions_memento['results']
        metabolights_dict['flags']['hasReactions'] = 'true'
        return metabolights_dict


class _FileHandler:

    @staticmethod
    def save_spectra(spectra_id, spectra_data, mtbls_id, destination) -> None:
        """
        Parse a given spectral data file into a dict, and then save it as a .json file using
        `FileHandler.save_json_file`
        :param spectra_id: Unique identifier of this spectral data file.
        :param spectra_data: spectral data file to process.
        :param mtbls_id: the MTBLC accession associated with this spectral data file.
        :param destination: The MTBLC directory to save the .json file to.
        :return: None
        """
        final_destination = f'{destination}/{mtbls_id}/{mtbls_id}_spectrum/{spectra_id}/{spectra_id}.json'
        datapoints = spectra_data.split(" ")
        ml_spectrum = {"spectrumId": spectra_id, "peaks": []}
        mz_array = []

        float_round_lambda = lambda num, places, direction: direction(num * (10 ** places)) / float(10 ** places)
        for datapoint in datapoints:
            temp_array = datapoint.split(":")
            temp_peak = {
                'intensity': float_round_lambda(float(temp_array[1].strip()) * 9.99, 6, math.floor),
                'mz': float(temp_array[0].strip())
            }
            ml_spectrum['peaks'].append(temp_peak)
            mz_array.append(float(temp_peak['mz']))
        ml_spectrum.update({'mzStart': min(mz_array), 'mzStop': max(mz_array)})
        _FileHandler.save_json_file(final_destination, ml_spectrum)

    @staticmethod
    def save_json_file(filename: str, data: dict) -> None:
        """
        Dump a given dict as a .json file. Check first that the directory we want to save to exists, and if it doesn't,
        create it.
        :param filename: string representation of the full path of the .json file to be.
        :param data: dict to be saved as a .json file
        :return: None
        """
        print(f'Attempting to save {filename}')
        if not os.path.exists(os.path.dirname(filename)):
            try:
                os.makedirs(os.path.dirname(filename))
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
        with open(filename, "w") as fp:
            try:
                json.dump(data, fp)
            except json.decoder.JSONDecodeError as e:
                print('what the hell ' + str(e))
        if os.path.exists(filename):
            print(f'Successfully saved {filename}')
        else:
            print(f'Failed to save {filename}')


class _InternalUtils:
    """
    Internal utils class, private to this script, sticking all static methods that don't belong anywhere else here (as
    util classes are SUPPOSED TO BE)
    """

    @staticmethod
    def get_val(root, key) -> str:
        """
        Just a try except wrapper method around an ET .find call. We do something similar with the xml_exception_angel
        annotation, but it seemed wrong to add that annotation to the entrypoint method for the chebi building portion
        of the script. Sue me.

        :param root: xml.etree.ElementTree object to search.
        :param key: The key to search the root document for.
        :return: Result of search
        """
        val = None
        try:
            val = root.find("{https://www.ebi.ac.uk/webservices/chebi}" + key).text
        except Exception as e:
            logging.exception(str(e))
            print(f'error getting key {str(e)} from chebi response')
        return val

    @staticmethod
    def initialize_compound_dict() -> dict:
        """
        Initialize the compound dict that will ultimately be saved. Pulled this almost wholesale from the old script,
        and chucked it into this method.

        :return: Initialized compound dict.
        """
        metabolights_compound = {}
        metabolights_compound["flags"] = {}
        metabolights_compound["flags"]['hasLiterature'] = "false"
        metabolights_compound["flags"]['hasReactions'] = "false"
        metabolights_compound["flags"]['hasSpecies'] = "false"
        metabolights_compound["flags"]['hasPathways'] = "false"
        metabolights_compound["flags"]['hasNMR'] = "false"
        metabolights_compound["flags"]['hasMS'] = "false"
        return metabolights_compound

    @staticmethod
    def preliminary_log_lines(metabolights_id) -> None:
        """
        Spit out the initial log lines.

        :param metabolights_id: The current accession ID.
        :return: None
        """
        _InternalUtils.gimme_line()
        print("Compound ID: " + metabolights_id)
        _InternalUtils.gimme_line()
        print("Process started: " + metabolights_id)
        print("Requesting compound chemical information from ChEBI:")

    @staticmethod
    def gimme_line():
        print('-----------------------------------------------------------------------------------------------------')

    @staticmethod
    def pascal_case(string: str) -> str:
        """
        Return a string with the first character lower case-d.

        :param string: String to pascal case-ify.
        :return: Pascal case-ified string.
        """
        return string[0].lower() + string[1:]

    @staticmethod
    def flag_is_enabled(rt_config: RuntimeFlags, method) -> bool:
        """
        Check whether a given config flag is enabled.
        :param rt_config: RuntimeFlags config object.
        :param method: Function wrapper that we extract the runtime flag from.
        :return: boolean indicating whether flag is enabled or not.
        """
        key = _InternalUtils.extract_name_from_function(method)
        string = rt_config.mapping[key]
        boolean = rt_config.__getattribute__(string)
        return boolean

    @staticmethod
    def extract_name_from_function(method) -> str:
        """
        Extract the exact function name from the function wrapper object.
        :param method: Function wrapper object.
        :return: Function name.
        """
        return method.__name__.split('.')[-1]
