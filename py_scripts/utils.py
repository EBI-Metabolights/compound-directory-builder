#!/usr/bin/env python
import errno
import logging
import os
import json
import threading
import time
#from urllib import request
import xml.etree.ElementTree as ET

#import urllib3.exceptions
from bs4 import BeautifulSoup
from math import ceil, floor
import datetime
import requests
import pickle

# MetaboLights api
from requests import Session

MetaboLightsWSUrl = "http://www.ebi.ac.uk/metabolights/ws/"
MetaboLightsWSStudyUrl = MetaboLightsWSUrl + "studies/public/study/"
MetaboLightsWSStudiesList = MetaboLightsWSUrl + "studies"
MetaboLightsWSCompoundsUrl = MetaboLightsWSUrl + "compounds/"
MetaboLightsWSCompoundsList = MetaboLightsWSUrl + "compounds/list"

# CHEBI api
chebiapi ="https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId="
chebiNSMap = {"envelop": "http://schemas.xmlsoap.org/soap/envelope/","chebi":"{http://www.ebi.ac.uk/webservices/chebi}"}

# Chemical Translation Service [ http://cts.fiehnlab.ucdavis.edu/ ]
ctsapi = "http://cts.fiehnlab.ucdavis.edu/service/compound/"

# CACTUS Chemical identifier resolver
cactusapi = "https://cactus.nci.nih.gov/chemical/structure/"

epmcAPI = "http://www.ebi.ac.uk/europepmc/webservices/rest/search?query="

MLStudiesList = []
MLCompoundsList = []

ReactomeJSON = {}
ReactomeUrl = "http://www.reactome.org/download/current/ChEBI2Reactome.txt"

#Rhea API
rheaapi = "https://www.rhea-db.org/rhea/"
chebiCompound = {}

class Timer:
    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.elapsed_time = self.end_time - self.start_time
        print(exc_type)
        print('Time elapsed:' + self.elapsed_time.__str__() + ' seconds')

def fetchMetaboLightsStudiesList():
    #studiesList = ["MTBLS626","MTBLS241","MTBLS762","MTBLS533","MTBLS1437"]
    #return studiesList
    return requests.get(MetaboLightsWSStudiesList).json()['content']

def fetchMetaboLightsCompoundsList():
    compound_list = []
    try:
        compound_list = requests.get(MetaboLightsWSCompoundsList).json()['content']
    except TimeoutError as e:
        logging.exception('Request timed out fetching compound list')
    except Exception as e:
        logging.exception('unexpected error')
    return compound_list

def generateMLStudyCompoundMappingFile(mappingFile, debug):
    logging.basicConfig(filename='mapping.log', level=logging.INFO)
    studiesList = {}
    compoundsList = {}
    speciesList = []
    mt = {}
    debug_counter = 0
    global MLStudiesList
    MLStudiesList = fetchMetaboLightsStudiesList()
    for study in MLStudiesList:
        debug_counter += 1
        if debug is True and debug_counter > 20:
            break
        print('processing ' + study)
        studyContent = requests.get(MetaboLightsWSStudyUrl + study).json()["content"]
        organismData = studyContent["organism"]
        hasMultipleOrganisms = True
        if organismData is None:
            hasMultipleOrganisms = False
        elif (len(organismData) == 1):
            hasMultipleOrganisms = False
        assayNumber = 1
        if studyContent["assays"] is None:
            continue
        for assay in studyContent["assays"]:
            try:
                mafFilePath = assay["metaboliteAssignment"]["metaboliteAssignmentFileName"]
                mafSplits = mafFilePath.split("/")
                mafFileName = mafSplits[len(mafSplits)-1]
                mafAPIurl = MetaboLightsWSStudiesList + "/" +study+ "/" +mafFileName
                metabolitesLines = requests.get(mafAPIurl).json()["data"]["rows"]
                #metabolitesLines = requests.get( MetaboLightsWSStudyUrl + study + "/assay/" + str(assayNumber) + "/maf").json()["content"]['metaboliteAssignmentLines']
                if metabolitesLines is None:
                    continue
                for line in metabolitesLines:
                    dbID = str(line['database_identifier'])
                    part = ""
                    if dbID != '':
                        species = str(line['species'])
                        if species == "":
                            if not hasMultipleOrganisms:
                                species = organismData[0]['organismName']
                                part = organismData[0]['organismPart']
                        if species not in speciesList and species != "":
                            speciesList.append(species)
                        tempCompound = {}
                        if dbID not in compoundsList:
                            compoundsList[dbID] = []
                            tempCompound['study'] = study
                            tempCompound['assay'] = assayNumber
                            tempCompound['species'] = species
                            tempCompound['part'] = part
                            tempCompound['taxid'] = line['taxid']
                            tempCompound['mafEntry'] = line
                            compoundsList[dbID].append(tempCompound)
                        else:
                            tempCompound['study'] = study
                            tempCompound['assay'] = assayNumber
                            tempCompound['species'] = species
                            tempCompound['part'] = part
                            tempCompound['taxid'] = line['taxid']
                            tempCompound['mafEntry'] = line
                            compoundsList[dbID].append(tempCompound)
                        tempStudy = {}
                        if study not in studiesList:
                            studiesList[study] = []
                            tempStudy['compound'] = dbID
                            tempStudy['assay'] = assayNumber
                            tempStudy['species'] = species
                            tempStudy['part'] = part
                            tempCompound['taxid'] = line['taxid']
                            studiesList[study].append(tempStudy)
                        else:
                            tempStudy['compound'] = dbID
                            tempStudy['assay'] = assayNumber
                            tempStudy['species'] = species
                            tempStudy['part'] = part
                            tempCompound['taxid'] = line['taxid']
                            studiesList[study].append(tempStudy)
                assayNumber += 1
            except Exception as ex:
                logging.exception('study ' + study + ' has been skipped')
                pass
    mt['studyMapping'] = studiesList
    mt['compoundMapping'] = compoundsList
    mt['speciesList'] = speciesList
    mt['updatedAt'] = str(datetime.datetime.now())
    print('attempting to save mapping.json file to ' + mappingFile)
    try:
        with Timer():
            with open(mappingFile, 'w') as mapping:
                json.dump(mt, mapping)
    except Exception as ex:
        print(str(ex))
        logging.exception('error in saving mapping.json file')

    if False:
        with Timer():
            save_via_pickle('mapping', mt)

def getReactomeData(reactomeFile):
    reactomeData = requests.get(ReactomeUrl).content
    for line in reactomeData.split("\n"):
        if line:
            dataArray = line.split("\t")
            metabolightsId = "MTBLC" + str(dataArray[0])
            tempDic = {}
            tempDic['reactomeId'] = str(dataArray[1])
            tempDic['reactomeUrl'] = str(dataArray[2])
            tempDic['pathway'] = str(dataArray[3])
            tempDic['pathwayId'] = str(dataArray[4])
            tempDic['species'] = str(dataArray[5])
        if metabolightsId not in ReactomeJSON:
            ReactomeJSON[metabolightsId] = [tempDic]
        else:
            ReactomeJSON[metabolightsId].append(tempDic)



    writeDataToFile( reactomeFile, ReactomeJSON)

def save_via_pickle(filename, obj):
    """Use when the json file is huge"""
    with open(f'{filename}.pickle', 'wb') as f:
        pickle.dump(obj, f)


def init(loggingObject):
	global logger
	logger = loggingObject

def getChebiData(chebiId,mlMapping):
    global chebiCompound
    print(chebiapi + chebiId)
    chebiRESTResponse = requests.get(chebiapi + chebiId).content
    root = ET.fromstring(chebiRESTResponse).find("envelop:Body", namespaces=chebiNSMap).find("{https://www.ebi.ac.uk/webservices/chebi}getCompleteEntityResponse").find("{https://www.ebi.ac.uk/webservices/chebi}return")
    print(root.find("{https://www.ebi.ac.uk/webservices/chebi}chebiId").text)
    # print root
    try:
        chebiCompound["id"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}chebiId").text
    except:
        pass
    
    try:
        chebiCompound["definition"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}definition").text
    except:
        pass

    try:
        chebiCompound["smiles"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}smiles").text
    except:
        pass
    try:
        chebiCompound["inchi"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}inchi").text
    except:
        pass

    try:
        chebiCompound["inchiKey"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}inchiKey").text
    except:
        pass

    try:
        chebiCompound["charge"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}charge").text
    except:
        pass
    
    try:
        chebiCompound["mass"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}mass").text
    except:
        pass

    try:
        chebiCompound["monoisotopicMass"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}monoisotopicMass").text
    except:
        pass

    try:
        chebiCompound["chebiAsciiName"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}chebiAsciiName").text
    except:
        pass

    try:
        chebiCompound["Synonyms"] = []
        for synonymn in root.findall("{https://www.ebi.ac.uk/webservices/chebi}Synonyms"):
            chebiCompound["Synonyms"].append(synonymn.find("{https://www.ebi.ac.uk/webservices/chebi}data").text)
    except:
        pass

    try:
        chebiCompound["IupacNames"] = []
        for iupacname in root.findall("{https://www.ebi.ac.uk/webservices/chebi}IupacNames"):
            chebiCompound["IupacNames"].append(synonymn.find("{https://www.ebi.ac.uk/webservices/chebi}data").text)
    except:
        pass

    try:
        chebiCompound["Formulae"] = root.find("{https://www.ebi.ac.uk/webservices/chebi}Formulae").find("{https://www.ebi.ac.uk/webservices/chebi}data").text
    except:
        pass

    try:
        chebiCompound["Citations"] = []
        for citation in root.findall("{https://www.ebi.ac.uk/webservices/chebi}Citations"):
            citationDic = {}
            citationDic["source"] = citation.find("{https://www.ebi.ac.uk/webservices/chebi}source").text
            citationDic["type"] = citation.find("{https://www.ebi.ac.uk/webservices/chebi}type").text
            citationDic["value"] = citation.find("{https://www.ebi.ac.uk/webservices/chebi}data").text
            chebiCompound["Citations"].append(citationDic)
    except:
        pass

    try:
        chebiCompound["DatabaseLinks"] = []
        for databaselink in root.findall("{https://www.ebi.ac.uk/webservices/chebi}DatabaseLinks"):
            databaselinkDic = {}
            databaselinkDic["source"] = databaselink.find("{https://www.ebi.ac.uk/webservices/chebi}type").text
            databaselinkDic["value"] = databaselink.find("{https://www.ebi.ac.uk/webservices/chebi}data").text
            chebiCompound["DatabaseLinks"].append(databaselinkDic)
    except:
        pass

    try:
        chebiCompound["CompoundOrigins"] = []
        chebiCompound["Species"] = {}
        for origin in root.findall("{https://www.ebi.ac.uk/webservices/chebi}CompoundOrigins"):
            chebispecies = origin.find("{https://www.ebi.ac.uk/webservices/chebi}speciesText").text.lower()
            if chebispecies not in chebiCompound["Species"]:
                originDic = {}
                chebiCompound["Species"][chebispecies] = []
                originDic["SpeciesAccession"] = origin.find("{https://www.ebi.ac.uk/webservices/chebi}speciesAccession").text
                originDic["SourceType"] = origin.find("{https://www.ebi.ac.uk/webservices/chebi}SourceType").text
                originDic["SourceAccession"] = origin.find("{https://www.ebi.ac.uk/webservices/chebi}SourceAccession").text
                chebiCompound["Species"][chebispecies].append(originDic)
            else:
                originDic = {}
                originDic["SpeciesAccession"] = origin.find("{https://www.ebi.ac.uk/webservices/chebi}speciesAccession").text
                originDic["SourceType"] = origin.find("{https://www.ebi.ac.uk/webservices/chebi}SourceType").text
                originDic["SourceAccession"] = origin.find("{https://www.ebi.ac.uk/webservices/chebi}SourceAccession").text
                chebiCompound["Species"][chebispecies].append(originDic)
    except:
        pass
    try:
        if chebiCompound["id"] in mlMapping['compoundMapping']:
            studyspecies = mlMapping['compoundMapping'][chebiCompound["id"]]
            for studyS in studyspecies:
                if (studyS['species'] != ""):
                    tempSSpecies = str(studyS['species']).lower()
                    if tempSSpecies not in chebiCompound["Species"]:
                        originDic = {}
                        chebiCompound["Species"][tempSSpecies] = []
                        originDic["Species"] = tempSSpecies
                        originDic["SpeciesAccession"] = studyS['study']
                        originDic["MAFEntry"] = studyS['mafEntry']
                        originDic["Assay"] = studyS['assay']
                        chebiCompound["Species"][tempSSpecies].append(originDic)
                    else:
                        originDic = {}
                        originDic["Species"] = tempSSpecies
                        originDic["SpeciesAccession"] = studyS['study']
                        originDic["MAFEntry"] = studyS['mafEntry']
                        originDic["Assay"] = studyS['assay']
                        chebiCompound["Species"][tempSSpecies].append(originDic)
    except:
        pass

def fetchCompound(metabolightsID, wd, dd, reactomeData, mlMapping, session: Session, file_lock: threading.Lock):
    logger.info("-----------------------------------------------")
    logger.info("Compound ID: " + metabolightsID)
    logger.info("-----------------------------------------------")
    logger.info("Process started: " + metabolightsID)
    logger.info("Requesting compound chemical information from ChEBI:")
    chebiId = metabolightsID.replace("MTBLC","").strip()
    mtblcs = session.get(MetaboLightsWSCompoundsUrl+metabolightsID).json()['content']
    #mtblcs = json.loads(request.urlopen(MetaboLightsWSCompoundsUrl+metabolightsID).read())["content"]
    getChebiData(chebiId, mlMapping)
    logger.info("Initialising MetaboLightsCompoundJSON")
    MetaboLightsCompoundJSON = {}
    MetaboLightsCompoundJSON["flags"] = {}
    MetaboLightsCompoundJSON["flags"]['hasLiterature'] = "false"
    MetaboLightsCompoundJSON["flags"]['hasReactions'] = "false"
    MetaboLightsCompoundJSON["flags"]['hasSpecies'] = "false"
    MetaboLightsCompoundJSON["flags"]['hasPathways'] = "false"
    MetaboLightsCompoundJSON["flags"]['hasNMR'] = "false"
    MetaboLightsCompoundJSON["flags"]['hasMS'] = "false"
    logger.info("Requesting compound chemical information from CTS:")
    try:
        if(chebiCompound["inchiKey"] != None):
            ctsc = session.get(ctsapi+chebiCompound["inchiKey"]).json()
            #ctsc = json.loads(request.urlopen(ctsapi+chebiCompound["inchiKey"]).read())
    except:
        pass
    MetaboLightsCompoundJSON["id"] = metabolightsID
    try:
        MetaboLightsCompoundJSON["name"] = chebiCompound["chebiAsciiName"]
    except:
        MetaboLightsCompoundJSON["name"] = "NA"
        logger.info("Compound Error: "+ metabolightsID + "Name not assigned")
        pass

    try:   
        MetaboLightsCompoundJSON["definition"] = chebiCompound["definition"]
    except:
        MetaboLightsCompoundJSON["definition"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Definition not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["iupacNames"] = chebiCompound["IupacNames"]
    except:
        MetaboLightsCompoundJSON["iupacNames"] = []
        logger.info("Compound Error: "+metabolightsID + "IUPAC Names not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["smiles"] = chebiCompound["smiles"]
    except:
        MetaboLightsCompoundJSON["smiles"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Smiles not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["inchi"] = chebiCompound["inchi"]
    except:
        MetaboLightsCompoundJSON["inchi"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Inchi not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["inchiKey"] = chebiCompound["inchiKey"]
    except:
        MetaboLightsCompoundJSON["inchiKey"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Inchikey not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["charge"] = chebiCompound["charge"]
    except:
        MetaboLightsCompoundJSON["charge"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Charge not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["averagemass"] = chebiCompound["mass"]
    except:
        MetaboLightsCompoundJSON["averagemass"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Average Mass not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["exactmass"] = chebiCompound["monoisotopicMass"]
    except:
        MetaboLightsCompoundJSON["exactmass"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Exact Mass not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["formula"] = chebiCompound["Formulae"]
    except:
        MetaboLightsCompoundJSON["formula"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Average Mass not assigned")
        pass

    logger.info("Fetching Citations:")
    try:
        MetaboLightsCompoundJSON["citations"] = getCitations(chebiCompound["Citations"], session)
        if MetaboLightsCompoundJSON["citations"]:
            MetaboLightsCompoundJSON["flags"]['hasLiterature'] = "true"
    except:
        MetaboLightsCompoundJSON["citations"] = []
        MetaboLightsCompoundJSON["flags"]['hasLiterature'] = "false"
        logger.info("Compound Error: "+metabolightsID + "Citations not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["species"] = chebiCompound["Species"]

        if MetaboLightsCompoundJSON["species"] :
            MetaboLightsCompoundJSON["flags"]['hasSpecies'] = "true"
    except:
        MetaboLightsCompoundJSON["species"] = []
        MetaboLightsCompoundJSON["flags"]['hasSpecies'] = "false"
        logger.info("Compound Error: "+metabolightsID + "Species not assigned")
        pass

    logger.info("Fetching Structure Data:")
    try:
        MetaboLightsCompoundJSON["structure"] = session.get(cactusapi+chebiCompound["inchiKey"]+"/sdf").text
        #MetaboLightsCompoundJSON["structure"] = request.urlopen(cactusapi+chebiCompound["inchiKey"]+"/sdf").read()
    except:
        MetaboLightsCompoundJSON["structure"] = "NA"
        logger.info("Compound Error: "+metabolightsID + "Structure not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["synonyms"] = chebiCompound["Synonyms"]
    except:
        MetaboLightsCompoundJSON["synonyms"] = []
        logger.info("Compound Error: "+metabolightsID + "Synonyms not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["pathways"] = mapPathways(
            chebiCompound, metabolightsID, reactomeData, MetaboLightsCompoundJSON, chebiCompound["inchiKey"], session)
    except:
        MetaboLightsCompoundJSON["pathways"] = {}
        logger.info("Compound Error: "+metabolightsID + "Pathways not assigned")
        pass

    try:
        MetaboLightsCompoundJSON["reactions"] = fetchReactions(chebiCompound, session)
        if MetaboLightsCompoundJSON["reactions"]:
            MetaboLightsCompoundJSON["flags"]['hasReactions'] = "true"
    except Exception as ex:
        logging.exception('there was some big error')
        MetaboLightsCompoundJSON["reactions"] = []
        MetaboLightsCompoundJSON["flags"]['hasReactions'] = "false"
        logger.info("Compound Error: "+ metabolightsID + "Reactions not assigned")
        pass
        
    MetaboLightsCompoundJSON["spectra"] = fetchSpectra(mtblcs['mc']['metSpectras'],metabolightsID, dd, session, logger)
    if MetaboLightsCompoundJSON["spectra"]['MS']:
        MetaboLightsCompoundJSON["flags"]['hasMS'] = "true"
    if MetaboLightsCompoundJSON["spectra"]['NMR']:
        MetaboLightsCompoundJSON["flags"]['hasNMR'] = "true"

    logger.info("Writing data to _data.json file - location: " + dd + "/" + metabolightsID + "/" + metabolightsID + "_data.json")
    file_lock.acquire()
    try:
        writeDataToFile(dd + "/" + metabolightsID + "/" + metabolightsID + "_data.json", MetaboLightsCompoundJSON)
    finally:
        file_lock.release()
    logger.info("-----------------------------------------------")

def fetchSpectra(spectra, metabolightsID, dd, session: Session, logger):
    MetSpec = {}
    MetSpec['NMR'] = []
    MetSpec['MS'] = []
    try:
        for spec in spectra:
            if(spec["spectraType"] == "NMR"):
                tempSpec = {}
                tempSpec["name"] = spec["name"]
                tempSpec["id"] = str(spec["id"])
                tempSpec["url"] = "http://www.ebi.ac.uk/metabolights/webservice/compounds/spectra/" + str(spec["id"]) + "/json"
                tempSpec["path"] = spec["pathToJsonSpectra"]
                tempSpec["type"] = spec["spectraType"]
                attriArray = []
                for attri in spec["attributes"]:
                    tempAttri = {}
                    tempAttri["attributeName"] = attri["attributeDefinition"]["name"]
                    tempAttri["attributeDescription"] = attri["attributeDefinition"]["description"]
                    tempAttri["attributeValue"] = attri["value"]
                    attriArray.append(tempAttri)
                tempSpec["attributes"] = attriArray
                MetSpec['NMR'].append(tempSpec)
    except:
        logger.info("Compound Error: "+ metabolightsID + " NMR Spectra not assigned")
        pass
    if "inchiKey" in chebiCompound.keys():
        MetSpec['MS'] = fetchMSFromMONA(chebiCompound["inchiKey"], metabolightsID, dd, session, logger)
   
    return MetSpec

def fetchMSFromMONA(inchikey, metabolightsID, dd, session: Session, logger):
    ml_spectrum = []
    url = "http://mona.fiehnlab.ucdavis.edu/rest/spectra/search?query=compound.metaData=q=%27name==\%22InChIKey\%22%20and%20value==\%22"+inchikey+"\%22%27"
    try:
        result = session.get(url).json()
    except requests.exceptions.RetryError as e:
        result = []
        logger.info('Unable to reach MoNA due to MaxRetryError: {0}'.format(str(e)))
    except requests.exceptions.ConnectionError as e:
        result = []
        logger.info('Unable to reach MoNA due to ConnectionError: {0}'.format(str(e)))
    except Exception as ex:
        result = []
        logging.exception('Error in trying to reach MoNA')
    #result = json.load(request.urlopen(url))
    for spectra in result:
        ml_spectra = {}
        ml_spectra['splash'] = spectra['splash']
        tempSpectraName = str(spectra['id'])
        ml_spectra['name'] =  tempSpectraName
        ml_spectra['type'] = "MS"
        ml_spectra['url'] = "/metabolights/webservice/beta/spectra/"+ metabolightsID + "/" + tempSpectraName
        tempSubmitter = spectra['submitter']
        ml_spectra['submitter'] =  str(tempSubmitter['firstName']) + " " + str(tempSubmitter['lastName']) + " ; " + str(tempSubmitter['emailAddress']) + " ; " + str(tempSubmitter['institution'])
        ml_spectra['attributes'] = []
        for metadata in spectra['metaData']:
            tempAttri = {}
            if metadata['computed'] == False:
                tempAttri['attributeName'] = metadata['name'] 
                tempAttri['attributeValue'] = metadata['value']
                tempAttri['attributeDescription'] = ""
                ml_spectra['attributes'].append(tempAttri)
        #if not copyright:
        ml_spectrum.append(ml_spectra)
        storeSpectra(tempSpectraName, spectra['spectrum'], metabolightsID, dd)
    return ml_spectrum

def storeSpectra(specid, spectraData, metabolightsID, dd):
    destination = dd + metabolightsID + "/" + metabolightsID + "_spectrum" + "/" + specid + "/" + specid + ".json"
    datapoints = spectraData.split(" ")
    mlSpectrum = {}
    mlSpectrum["spectrumId"] = specid
    mlSpectrum["peaks"] = []
    mzArray = []
    for datapoint in datapoints:
        tempArray = datapoint.split(":")
        tempPeak = {}
        tempPeak['intensity'] = float_round(float(tempArray[1].strip()) * 9.99, 6)
        tempPeak['mz'] = float(tempArray[0].strip())
        mlSpectrum["peaks"].append(tempPeak)
        mzArray.append(float(tempPeak['mz']))
    mlSpectrum["mzStart"] = min(mzArray)
    mlSpectrum["mzStop"] = max(mzArray)
    writeDataToFile(destination, mlSpectrum)

def float_round(num, places = 0, direction = floor):
    return direction(num * (10**places)) / float(10**places)

def fetchReactions(chebi, session: Session):
    query = "?query="
    columns = "&columns=rhea-id,equation,chebi-id"
    format = "&format=json"
    limit = "&limit=10"
    rheaData = session.get(rheaapi + query + chebi["id"] + columns + format + limit).json()
    if len(rheaData['results']) is 0:
        return []
    reactions = []

    for result in rheaData["results"]:
        reaction = {}
        reaction['name'] = result["equation"]
        reaction["id"] = result["id"]
        # two below fields used to be populated responses from API's that no longer exist but keeping them for posterity
        reaction["biopax2"] = ""
        reaction["cmlreact"] = ""
        reactions.append(reaction)

    return reactions

def mapPathways(chebi, compound, reactomeData, MetaboLightsCompoundJSON, InChIKey, session: Session):
    tempPathwayDictionary = {}
    tempPathwayDictionary["WikiPathways"] = getWikiPathwaysData(InChIKey, session)
    tempPathwayDictionary["KEGGPathways"] = getKEGGData(chebi, session)
    tempPathwayDictionary["ReactomePathways"] = getReactomePathwaysData(compound, reactomeData)
    
    if (len(tempPathwayDictionary["WikiPathways"]) > 0 ):
        MetaboLightsCompoundJSON["flags"]['hasPathways'] = "true"
    elif(len(tempPathwayDictionary["ReactomePathways"]) >0 ):
        MetaboLightsCompoundJSON["flags"]['hasPathways'] = "true"
    elif(len(tempPathwayDictionary["KEGGPathways"]) > 0):
        MetaboLightsCompoundJSON["flags"]['hasPathways'] = "true"

    return tempPathwayDictionary

def getWikiPathwaysData(InChIKey, session: Session):
    pathways = {}
    try:
        wikipathwaysapi = "https://webservice.wikipathways.org/findPathwaysByXref?ids="+InChIKey+"&codes=Ik&format=json"
        #wikipathwaysapi = "http://webservice.wikipathways.org/index.php?ids="+chebi["id"]+"&codes=Ce&method=findPathwaysByXref&format=json"
        wikipathways = session.get(wikipathwaysapi).json()['result']
        #wikipathways =  json.loads(request.urlopen(wikipathwaysapi).read())["result"]
        if len(wikipathways) > 0 :
            for pathway in wikipathways:
                if pathway["species"] in pathways:
                    tempDict = {}
                    tempDict["id"] = pathway["id"]
                    tempDict["url"] = pathway["url"]
                    tempDict["name"] = pathway["name"]
                    if tempDict not in pathways[pathway["species"]]:
                        pathways[pathway["species"]].append(tempDict)
                else:
                    pathways[pathway["species"]] = []
                    tempDict = {}
                    tempDict["id"] = pathway["id"]
                    tempDict["url"] = pathway["url"]
                    tempDict["name"] = pathway["name"]
                    pathways[pathway["species"]].append(tempDict)
    except:
        pass
    return pathways
    
def getReactomePathwaysData(compound, reactomeData):
    tempReactomePathways = reactomeData[compound]
    reactomePathways = {}
    try:
        for pathway in tempReactomePathways:
            tempPathway = {}
            tempPathway['name'] = pathway['pathway']
            tempPathway['pathwayId'] = pathway['pathwayId']
            tempPathway['url'] = pathway['reactomeUrl']
            tempPathway['reactomeId'] = pathway['reactomeId']
            if pathway['species'] not in reactomePathways:
                reactomePathways[pathway['species']] = [tempPathway]
            else:
                reactomePathways[pathway['species']].append(tempPathway)
    except:
        pass
    return reactomePathways

def getKEGGData(chebi, session: Session):
    pathways = []
    try:
        keggapi = "http://rest.kegg.jp/conv/compound/" + chebi["id"].lower()
        keggid = session.get(keggapi).text.split("\t")[1].strip()
        #keggid = request.urlopen(keggapi).read().split("\t")[1].strip()
        pathwayslistapi = "http://rest.kegg.jp/link/pathway/" + keggid
        pathwaysData = session.get(pathwayslistapi).text
        # pathwaysData = request.urlopen(pathwayslistapi).read()
        for line in pathwaysData.strip().split("\n"):
            pathwayId = line.split("\t")[1].strip()
            pathwayapi = "http://rest.kegg.jp/get/" + pathwayId
            pathwayData = session.get(pathwayapi).text
            #pathwayData = request.urlopen(pathwayapi).read()
            tempDict = {}
            tempDict["id"] = pathwayId
            for pline in pathwayData.strip().split("\n"):
                if "NAME" in pline:
                    tempDict["name"] = pline.replace("NAME","").strip()
                elif "KO_PATHWAY" in pline:
                    tempDict["KO_PATHWAY"] = pline.replace("KO_PATHWAY","").strip()
                elif "DESCRIPTION" in pline:
                    tempDict["description"] = pline.replace("DESCRIPTION","").strip()
            pathways.append(tempDict)
    except:
        pass
    return pathways

def mapExternalIDS(chebi,ctsc):
    externalIDs = {}
    for link in chebi["DatabaseLinks"]:
        if link["source"] != "":
            if link["source"] in externalIDs:
                externalIDs[link["source"]].append(link["value"])
            else:
                externalIDs[link["source"]] = [link["value"]]

    for link in ctsc["externalIds"]:
        if link["name"] != "":
            if link["name"] in externalIDs:
                externalIDs[link["name"]].append(link["value"])
            else:
                externalIDs[link["name"]] = [link["value"]]
    return externalIDs

def getSynonymns(chebi,ctsc):
    synonyms = chebi["Synonyms"]
    # # cts synonyms
    # for synonym in ctsc["synonyms"]:
    #     synonyms.append(synonym["name"])
    # # remove duplicates
    # synonyms = list(set(name.lower() for name in synonyms))
    return synonyms

def getCitations(citationsList, session: Session):
    epmcList = []
    for citation in citationsList:
        try:
            tempCitation = citation
            epmc = epmcAPI + str(citation['value']) + "&format=json&resulttype=core"
            epmcData = session.get(epmc)
            #epmcData = request.urlopen(epmc).read()
            citationData =  epmcData.json()['resultList']['result'][0]
            tempCitation['title'] = citationData['title']
            tempCitation['abstract'] = citationData['abstractText']
            try:
                tempCitation['doi'] = citationData['doi']
            except:
                tempCitation['doi'] = "NA"
                pass
            tempCitation['author'] = citationData['authorString']
            epmcList.append(tempCitation)
        except:
            pass
    return epmcList

def getReactomeData(reactomeFile):
    session = Session()
    resp = session.get(ReactomeUrl)
    #reactomeData = urllib.urlopen(ReactomeUrl).read()
    for line in resp.text.split("\n"):
        if line:
            dataArray = line.split("\t")
            metabolightsId = "MTBLC" + str(dataArray[0])
            tempDic = {}
            tempDic['reactomeId'] = str(dataArray[1])
            tempDic['reactomeUrl'] = str(dataArray[2])
            tempDic['pathway'] = str(dataArray[3])
            tempDic['pathwayId'] = str(dataArray[4])
            tempDic['species'] = str(dataArray[5])
        if metabolightsId not in ReactomeJSON:
            ReactomeJSON[metabolightsId] = [tempDic]
        else:
            ReactomeJSON[metabolightsId].append(tempDic)

    print('attempting to save reactone.json file to ' + reactomeFile)
    with open(reactomeFile, "w") as fp:
        json.dump(ReactomeJSON, fp)

def getDateAndTime():
	ts = time.time()
	return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

def writeDataToFile(filename, data):
	if not os.path.exists(os.path.dirname(filename)):
		try:
			os.makedirs(os.path.dirname(filename))
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
	with open(filename, "w") as fp:
		json.dump(data, fp)

def appendDataToFile(filename, data):
	if not os.path.exists(os.path.dirname(filename)):
		try:
			os.makedirs(os.path.dirname(filename))
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
	with open(filename, "a") as fp:
		json.dump(data, fp)
