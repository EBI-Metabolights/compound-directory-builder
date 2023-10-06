from pydantic import BaseModel


class RuntimeFlags(BaseModel):
    spectra: bool = False
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
    #new_mona_api: str = "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search?query=exists(compound.metaData.name%3A'InChIKey'%20and%20compound.metaData.value%3A'{0}')"
    new_mona_api: str = "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search?endpoint=search&query=exists(compound.metaData.name%3A%27InChIKey%27%20and%20compound.metaData.value%3A%27{0}%27)%20and%20exists((tags.text%3A%27LC-MS%27))&size=10"

class CompoundBuilderUrls(BaseModel):
    mtbls: MtblsWsUrls = MtblsWsUrls()
    kegg: KeggUrls = KeggUrls()
    misc_urls: MiscUrls = MiscUrls()


class CompoundBuilderConfig:
    objs: CompoundBuilderObjs = CompoundBuilderObjs()
    urls: CompoundBuilderUrls = CompoundBuilderUrls()
    rt_flags: RuntimeFlags = RuntimeFlags()

