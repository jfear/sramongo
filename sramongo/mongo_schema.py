"""Set up MonoDB schema using monogodngine."""
from textwrap import fill
import datetime

from mongoengine import Document, EmbeddedDocument
from mongoengine import StringField, IntField, FloatField, \
    ListField, DictField, MapField, DateTimeField
from mongoengine import EmbeddedDocumentField, signals
from mongoengine.errors import ValidationError, FieldDoesNotExist

from sramongo.sra import SraExperiment
from sramongo.logger import logger


def handler(event):
    """Signal decorator to allow use of callback functions."""
    def decorator(fn):
        def apply(cls):
            event.connect(fn, sender=cls)
            return cls

        fn.apply = apply
        return fn
    return decorator


@handler(signals.pre_save)
def update_modified(sender, document):
    """Update db_modified before saving."""
    document.db_modified = datetime.datetime.now()


class DocumentString(object):
    def __init__(self, document):
        """Represents mongo Document as string.

        This is a helper class to make a string representation of a mongo
        document. In other words allows pretty printing of documents.

        Parameters
        ----------
        document: mongoengine.Document
            A subclass of the mongoengine.Document. Assumes that all data is
            stored in an attribute self._data.
        """
        self.string = ''
        self.add_dict(document._data)

    def add_dict(self, dictionary, spacer=0):
        """Dictionary string representation.

        Iterates over a dictionary and builds key: value string
        representations. Carries indentation level through other methods so that
        hierarchy is visually maintained.

        Parameters
        ----------
        dictionary: dict
            A dictionary with values to be turned into a string.

        spacer: int
            The amount of indention to use. As we go down further into a
            hierarchy indent values further.
        """
        for k, v in dictionary.items():
            if k.startswith('_'):       # skip private keys
                continue
            if isinstance(v, dict):     # if another dict call itself
                self.string += ' ' * spacer + '{}:\n'.format(k)
                self.add_dict(v, spacer=spacer + 4)
            elif isinstance(v, str) | (v is None):
                self.add_element(k, v, spacer=spacer)
            elif isinstance(v, list):
                self.add_list(k, v, spacer=spacer)
            elif isinstance(v, EmbeddedDocument) | isinstance(v, Document):
                self.add_embeded_document(k, v, spacer=spacer)
            else:
                print(k, v, type(v))
                pass

    def add_element(self, key, value, spacer=0):
        """Adds key: value pair when value is a string.

        Parameters
        ----------
        key: str
            Key of the item.

        value: str
            Value of the item.

        spacer: int
            Length of the indention.

        """
        spacer = ' ' * spacer
        if value is not None:
            # wrap really long lines like abstracts
            if len(value) > 100:
                s = spacer + ' ' * (len(key) + 2)
                value = fill(value, width=75).replace('\n', '\n' + s)
        self.string += spacer + '{}: {}\n'.format(key, value)

    def add_list(self, key, value, spacer=0):
        """Adds key: value pair when value is a list.

        Parameters
        ----------
        key: str
            Key of the item.

        value: list
            Value of the item.

        spacer: int
            Length of the indention.

        """
        spacer = ' ' * spacer
        if len(value) > 0:
            self.string += spacer + '{}:\n'.format(key)
            for x in value:
                s = spacer + ' ' * 4
                v = str(x).strip().replace('\n', ', ')
                self.string += s + '{}\n'.format(v)
        else:
            self.string += spacer + '{}: None\n'.format(key)

    def add_embeded_document(self, key, value, spacer=0):
        """Adds key: value pair when value is an embedded document.

        Parameters
        ----------
        key: str
            Key of the item.

        value: mongoengine.Document|mongoengine.EmbeddedDocument
            Value of the item.

        spacer: int
            Length of the indention.

        """
        spacer = ' ' * spacer
        self.string += spacer + '{}:\n'.format(key)
        s = spacer + ' ' * 4
        values = str(value).strip().split('\n')
        v = '\n'.join([s + x for x in values]) + '\n'
        self.string += v


# Different Types of Links
class URLLink(EmbeddedDocument):
    label = StringField()
    url = StringField()


class Xref(EmbeddedDocument):
    db = StringField()
    id = StringField()


class XrefLink(EmbeddedDocument):
    label = StringField()
    db = StringField()
    id = StringField()
    meta = {'allow_inheritance': True}


class EntrezLink(XrefLink):
    query = StringField()


class DDBJLink(XrefLink):
    url = StringField()


class ENALink(XrefLink):
    url = StringField()


# Organization
class Organization(EmbeddedDocument):
    """Organization embedded document.

    An organization contains information about the group that submitted to sra.
    For example, all data submitted to GEO are submitted to SRA using the GEO
    credentials.

    Attributes
    ----------
    organization_type: str
        Weather this organization is a center or individual or some other kind
        of group.

    abbreviation: str
        A short name for the organization.

    name: str
        Name of the organization.

    emai: str
        Contact email address.

    first_name: str
        First name of the person who submitted the data.

    last_name: str
        First name of the person who submitted the data.

    """
    organization_type = StringField()
    abbreviation = StringField()
    name = StringField()
    email = StringField()
    first_name = StringField()
    last_name = StringField()


# Submission
class Submission(EmbeddedDocument):
    """Submission embedded document.

    A submission is a group of experiments that are uploaded. This is the
    highest level container for information in the SRA. Typically not useful
    unless trying to track down the submitter.

    Attributes
    ----------
    submission_id: mongoengine.StringField
        This is ID will start with SRA/ERA/DRA depending on the originating
        database.

    broker: mongoengine.StringField
        Name of the group that submitted the data.

    external_id: mongoengine.ListField of mongoengine.EmbeddedDocumentField(Xref)
        External identifiers stored as {'db': value, 'id': value}

    secondary_id: mongoengine.ListField of mongoengine.EmbeddedDocumentField(Xref)
        Secondary identifiers stored as {'db': value, 'id': value}

    submitter_id: mongoengine.ListField of mongoengine.EmbeddedDocumentField(Xref)
        submitter identifiers stored as {'db': value, 'id': value}

    uuid: mongoengine.ListField
        A list of associated uuids.

    """
    submission_id = StringField()
    broker = StringField()

    # Other IDS
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)
    uuid = ListField(StringField(), default=list)


# Study
class RelatedStudy(XrefLink):
    is_primary = StringField()


class Study(EmbeddedDocument):
    """The contents of a SRA study.

    A study consists of a set of experiments designed with an overall goal in
    mind. For example, this could include a control experiment and a treatment
    experiment with the goal being to identify expression differences resulting
    from the treatment. The SRA study is the top level of the submission
    hierarchy.

    Attributes
    ----------
    study_id: mongoengine.StringField
        The primary identifier for a study. Identifiers begin with
        SRP/ERP/DRP depending on which database they originate from.

    GEO: mongoengine.StringField
        An identifier relating a study to the GEO database.

    GEO_Dataset: mongoengine.StringField
        An identifier relating to a study to the corresponding GEO datasets
        database.

    BioProject: mongoengine.StringField
        An identifier relating a study to the BioProject.

    pubmed: mongoengine.StringField
        An identifier relating a study to any published papers that used this
        study.

    external_id: mongoengine.ListField
        List of additional external ids.

    secondary_id: mongoengine.ListField
        List of secondary ids.

    submitter_id: mongoengine.ListField
        List of submitter ids.

    uuid: mongoengine.ListField
        List of uuids.

    title: mongoengine.StringField
        The title of the study.

    study_type: mongoengine.StringField
        The value is one of :ref:`studyTypes`

    abstract: mongoengine.StringField
        Abstract of the study.

    center_name: mongoengine.StringField
        Name of the submitting center.

    center_project_name: mongoengine.StringField
        Center specific identifier for the study.

    description: mongoengine.StringField
        Additional text describing the study.

    related_studies: mongoengine.ListField
        List of related studies.

    url_links: mongoengine.ListField
        List of url links

    xref_links: mongoengine.ListField
        List of external database cross references.

    entrez_links: mongoengine.ListField
        List of NCBI Entrez links.

    ddbj_links: mongoengine.ListField
        List of Japan's ddbj links.

    ena_links: mongoengine.ListField
        List of Europe's ena links.

    """
    study_id = StringField()

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    GEO_Dataset = StringField()
    BioProject = StringField()
    pubmed = StringField()

    # Other IDs
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)
    uuid = ListField(StringField(), default=list)

    # Attributes
    title = StringField()
    study_type = StringField()
    abstract = StringField()
    center_name = StringField()
    center_project_name = StringField()
    description = StringField()

    # Related study info
    related_studies = ListField(EmbeddedDocumentField(RelatedStudy), default=list)

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)


class Attribute(EmbeddedDocument):
    name = StringField()
    value = StringField()


# Sample
class Sample(EmbeddedDocument):
    """The contents of a SRA sample.

    A sample is the biological unit. An individual sample or a pool of samples
    can be use in the SRA Experiment. This document contains information
    describing the sample ranging from species information to detailed
    descriptions of what and how material was collected.

    Attributes
    ----------
    sample_id: mongoengine.StringField
        The primary identifier for a sample. Identifiers begin with
        SRS/ERS/DRS depending on which database they originate from.

    GEO: mongoengine.StringField
        An identifier relating a sample to the GEO database.

    BioSample: mongoengine.ReferenceField
        References the corresponding BioSample.

    BioProject: mongoengine.StringField
        An identifier relating a sample to the BioProject.

    pubmed: mongoengine.StringField
        An identifier relating a study to any published papers that use this
        sample.

    external_id: mongoengine.ListField
        List of additional external ids.

    secondary_id: mongoengine.ListField
        List of secondary ids.

    submitter_id: mongoengine.ListField
        List of submitter ids.

    uuid: mongoengine.ListField
        List of uuids.

    title: mongoengine.StringField
        The title of the sample.

    taxon_id: mongoengine.StringField
        The NCBI taxon id.

    scientific_name: mongoengine.StringField
        The scientific name.

    common_name: mongoengine.StringField
        The common name.

    individual_name: mongoengine.StringField
        The sample name.

    description: mongoengine.StringField
        Additional text describing the sample.

    attributes: mongoengine.DictField
        A set of key:value pairs describing the sample. For example tissue:ovary
        or sex:female.

    url_links: mongoengine.ListField
        List of url links

    xref_links: mongoengine.ListField
        List of external database cross references.

    entrez_links: mongoengine.ListField
        List of NCBI Entrez links.

    ddbj_links: mongoengine.ListField
        List of Japan's ddbj links.

    ena_links: mongoengine.ListField
        List of Europe's ena links.

    """
    # SRS/DRS/ERS
    sample_id = StringField()

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    BioSample = StringField()
    BioProject = StringField()
    pubmed = StringField()

    # Other IDs
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)
    uuid = ListField(StringField(), default=list)

    # Attributes
    title = StringField()
    taxon_id = StringField()
    scientific_name = StringField()
    common_name = StringField()
    individual_name = StringField()
    description = StringField()
    attributes = ListField(EmbeddedDocumentField(Attribute), default=list)

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)


# Experiment
class Experiment(EmbeddedDocument):
    """Experiment Document.

    An experiment describes an individual library. This library is made from a
    biological sample. Multiple experiments (or libraries) make up a study. An
    experiment can be sequenced on multiple lanes making up different Runs.

    Attributes
    ----------
    experiment_id: mongoengine.StringField
        The primary identifier for an experiment. Identifiers begin with
        SRX/ERX/DRX depending on which database they originate from.

    GEO: mongoengine.StringField
        An identifier relating an experiment to the GEO database.

    GEO_Dataset: mongoengine.StringField
        An identifier relating to an experiment to the corresponding GEO
        datasets database.

    pubmed: mongoengine.StringField
        An identifier relating an experiment to any published papers.

    external_id: mongoengine.ListField
        List of additional external ids.

    secondary_id: mongoengine.ListField
        List of secondary ids.

    submitter_id: mongoengine.ListField
        List of submitter ids.

    uuid: mongoengine.ListField
        List of uuids.

    title: mongoengine.StringField
        The title of the experiment.

    study_id: mongoengine.StringField
        The id for the corresponding study.

    design: mongoengine.StringField
        Free text description of the design of the experiment. Really anything
        could be in this field.

    library_name: mongoengine.StringField
        A identifier for this experiment.

    library_strategy: mongoengine.StringField
        The type of library constructed, one of :ref:`libraryStrategy`.

    library_source: mongoengine.StringField
        The source of library, one of :ref:`librarySource`.

    library_selection: mongoengine.StringField
        The selection method of library, one of :ref:`librarySelection`.

    library_layout: mongoengine.StringField
        The layout of the library, paired or single ended.
        :ref:`libraryLayout`.

    library_layout_orientation: mongoengine.StringField

    library_layout_length: mongoengine.StringField
        The average length of reads in this library.

    library_layout_sdev: mongoengine.StringField
        The standard deviation of length of reads in this library.

    pooling_strategey: mongoengine.StringField
        A description of how samples were pooled.

    library_construction_protocol: mongoengine.StringField
        A free text description of how the libraries were constructed.

    platform: mongoengine.StringField
        The platform used for sequencing, one of :ref:`platforms`.

    instrument_model: mongoengine.StringField
        The instrument model, one of :ref:`instrumentModels`.

    attributes: mongoengine.ListField of mongoengine.DictField
        A list of dictionaries containing key:value pairs describing the
        experiment. The stored dictionaries are of the form {'name': value,
        'value': value}. This was done to make querying easier.

    url_links: mongoengine.ListField
        List of url links

    xref_links: mongoengine.ListField
        List of external database cross references.

    entrez_links: mongoengine.ListField
        List of NCBI Entrez links.

    ddbj_links: mongoengine.ListField
        List of Japan's ddbj links.

    ena_links: mongoengine.ListField
        List of Europe's ena links.

    """
    # SRX/DRX/ERX
    experiment_id = StringField()

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    GEO_Dataset = StringField()
    pubmed = StringField()

    # Other IDs
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)
    uuid = ListField(StringField(), default=list)

    # Attributes
    title = StringField()
    study_id = StringField()
    design = StringField()
    library_name = StringField()
    library_strategy = StringField()
    library_source = StringField()
    library_selection = StringField()
    library_layout = StringField()
    library_layout_orientation = StringField()
    library_layout_length = StringField()
    library_layout_sdev = StringField()
    pooling_strategey = StringField()
    library_construction_protocol = StringField()
    platform = StringField()
    instrument_model = StringField()
    attributes = ListField(EmbeddedDocumentField(Attribute), default=list)

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)


# Run
class TaxRecord(EmbeddedDocument):
    parent = StringField()
    total_count = IntField()
    self_count = IntField()
    tax_id = StringField()
    name = StringField()


class TaxAnalysis(EmbeddedDocument):
    nspot_analyze = IntField()
    total_spots = IntField()
    mapped_spots = IntField()
    tax_counts = MapField(ListField(EmbeddedDocumentField(TaxRecord), default=list))


class Run(EmbeddedDocument):
    """Run Document.

    A Run describes a dataset generated from an Experiment. For example if a
    Experiment is sequenced on multiple lanes of a Illumina flowcell then data
    from each lane are considered a Run.

    Attributes
    ----------
    run_id: mongoengine.StringField
        The primary identifier for a run. Identifiers begin with
        SRR/ERR/DRR depending on which database they originate from.

    external_id: mongoengine.ListField
        List of additional external ids.

    secondary_id: mongoengine.ListField
        List of secondary ids.

    submitter_id: mongoengine.ListField
        List of submitter ids.

    uuid: mongoengine.ListField
        List of uuids.

    experiment_id: mongoengine.StringField
        The ID of the corresponding experiment.

    samples: mongoengine.ListField
        A list of sample IDs.

    nspots: mongoengine.IntField
        The total number of spots on a Illumina flowcell.

    nbases: mongoengine.IntField
        The number of bases.

    tax_analysis = EmbeddedDocumentField(TaxAnalysis)
        A dictionary containing results from a taxonomic analysis. Some Runs are
        analyzed and the number of reads that align to different taxa are
        recorded. The taxanomic analysis is stored in the SRA as a hierarchy,
        but it is stored here as a flattend dictionary for easier access to
        different classes. Basic structure is:

            'nspoot_analyze': The number of spots analyzed,
            'total_spots': The total number of spots,
            'mapped_spots': The number of spots that were able to be mapped,
            'tax_count': A dictionary containing actual taxonomic counts
            organized by level in the tree of life

                'kingdom':
                    ...
                'species':
                    'parent':
                        Name of parent level.
                    'total_count':
                        Number of mapped spots at this level and below.
                    'self_count':
                        Number of mapped spots at this level.
                    'tax_id':
                        taxonomic identifier.
                    'name':
                        of this taxonomy.
                'subspeciies':
                    ...

    nreads: mongoengine.IntField
        The number of reads.

    read_count_r1: mongoengine.FloatField
        Some Runs have additional information on reads. This is the number of
        reads from single ended or the first read pair in pair ended data.

    read_len_r1: mongoengine.FloatField
        This is the average length of reads from single ended or the first read
        pair in pair ended data.

    read_count_r2: mongoengine.FloatField
        This is the number of reads from the second read pair in pair ended
        data.

    read_len_r2: mongoengine.FloatField
        This is the avearge length of reads from the second read pair in pair
        ended data.

    release_date: mongoengine.StringField
        Release date of the Run. This information is from the runinfo table and
        not the XML.

    load_date: mongoengine.StringField
        Date the Run was uploaded. This information is from the runinfo table
        and not the XML.

    size_MB: mongoengine.IntField
        Size of the Run file. This information is from the runinfo table and not
        the XML.

    download_path: mongoengine.StringField
        Download path of the Run file. This information is from the runinfo
        table and not the XML.

    run_flags: mongoengine.ListField
        These are custom flags that I add. I have created these flags by
        interpreting specific information from data provided by the SRA. A list
        of :ref:`database_flags` and their descriptions.

    """
    # SRR/DRR/ERR
    run_id = StringField()

    # Other IDs
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)
    uuid = ListField(StringField(), default=list)

    # Attributes
    experiment_id = StringField()
    samples = ListField(StringField(), default=list)
    nspots = IntField()
    nbases = IntField()
    tax_analysis = EmbeddedDocumentField(TaxAnalysis)
    nreads = IntField()

    # if single ended then just use _r1
    read_count_r1 = FloatField()
    read_len_r1 = FloatField()

    read_count_r2 = FloatField()
    read_len_r2 = FloatField()

    # NOTE: This field is added by me to summarize
    run_flags = ListField(StringField(), default=list)

    # NOTE: Additional Fields not in the SRA XML but in summary table
    release_date = StringField()
    load_date = StringField()
    size_MB = IntField()
    download_path = StringField()


# SRA holder subdocument
class Sra(EmbeddedDocument):
    """SubDocument containg all SRA data.

    This is the general holder class for all SRA data. It has a couple of
    summary fields generate at import time.

    db_flags: mongoengine.ListField
        List of :ref:`database_flags`.

    db_imported: mongoengine.DateTimeField
        The date this document was added to the mongo database. In other words
        when you downloaded the data.

    """

    submission = EmbeddedDocumentField(Submission)
    organization = EmbeddedDocumentField(Organization)
    study = EmbeddedDocumentField(Study)
    sample = EmbeddedDocumentField(Sample)
    experiment = EmbeddedDocumentField(Experiment)
    run = ListField(EmbeddedDocumentField(Run), default=list)
    pool = ListField(StringField(), default=list)
    db_flags = ListField(StringField(), default=list)
    db_imported = DateTimeField(default=datetime.datetime.now)


# BioSample
class Contacts(EmbeddedDocument):
    email = StringField()
    first_name = StringField()
    last_name = StringField()


class BioSample(EmbeddedDocument):
    """The contents of a BioSample.

    BioSample is another database housed at NCBI which records sample metadata.
    This information should already be present in the Sra.sample information,
    but to be safe we can pull into the BioSample for additional metadata.

    Attributes
    ----------
    biosample_accn: mongoengine.StringField
        The primary identifier for a BioSample. These are the accession number
        which begin with SAM.

    biosample_id: mongoengine.StringField
        The primary identifier for a BioSample. These are the id number.

    biosample_primary: mongoengine.StringField
        This is the identifier that was said to be the primary identifier for a
        BioSample. These appear to be just the accession with SAM.

    biosample_secondary: mongoengine.StringField
        A secondary identifier for a BioSample.

    sample_id: mongoengine.StringField
        Unique id of the sample provided by the submitter.

    GEO: mongoengine.StringField
        GEO sample ID (GSM).

    title: mongoengine.StringField
        A free text description of the sample.

    description: mongoengine.StringField
        A free text description of the sample.

    tax_id: mongoengine.StringField
        The tax_id that the sample belongs.

    tax_name: mongoengine.StringField
        Scientific name of the organism that sample is from.

    organism_name: mongoengine.StringField
        Scientific name of the organism that sample is from.

    institute: mongoengine.StringField
        Name of the submitting institute.

    access: mongoengine.StringField
        Type of access, typically "public".

    publication_date: mongoengine.StringField
        Date the sample was published.

    last_update: mongoengine.StringField
        Last time BioSample updated sample information.

    submission_date: mongoengine.StringField
        Date the sample was submitted

    contacts: mongoengine.EmbeddedDocumentField
        Dictionary of contact information.

    models: mongoengine.ListField
        List of model information.

    attributes: mongoengine.ListField of mongoengine.DictField
        A list of dictionaries containing key:value pairs describing the
        experiment. The stored dictionaries are of the form {'name': value,
        'value': value}. This was done to make querying easier.

    """
    biosample_accn = StringField()
    biosample_id = StringField()
    biosample_primary = StringField()
    biosample_secondary = StringField()
    sample_id = StringField()
    GEO = StringField()
    title = StringField()
    description = StringField()
    tax_id = StringField()
    tax_name = StringField()
    organism_name = StringField()
    institute = StringField()
    access = StringField()
    publication_date = StringField()
    last_update = StringField()
    submission_date = StringField()
    contacts = ListField(EmbeddedDocumentField(Contacts), default=list)
    models = ListField(StringField())
    attributes = ListField(EmbeddedDocumentField(Attribute), default=list)


# BioProject
class BioProject(EmbeddedDocument):
    """The contents of a BioProject.

    BioProject is another database housed at NCBI which records project
    metadata.  This information should already be present in the SRA
    information, but to be safe we can pull into the BioProject for additional
    metadata.

    Attributes
    ----------
    bioproject_accn: mongoengine.StringField
        The primary identifier for a BioProject. These are the accession number
        which begin with PRJ.

    bioproject_id: mongoengine.StringField
        The primary identifier for a BioProject. These are the id numbers.

    name: mongoengine.StringField
        A brief name of the project.

    title: mongoengine.StringField
        The title of the project.

    description: mongoengine.StringField
        A short description of the project.

    publication: mongoengine.StringField
        Publication information.

    publication_date: mongoengine.StringField
        Date of publication.

    submission_id: mongoengine.StringField
        Identifier for the submission.

    last_date: mongoengine.DateTimeField
        Last date the BioProject was updated.

    submission_date: mongoengine.DateTimeField
        Date the BioProject was submitted.

    external_id: mongoengine.ListField
        List of additional external ids.

    """
    bioproject_accn = StringField()
    bioproject_id = StringField()
    name = StringField()
    title = StringField()
    description = StringField()
    publication = StringField()
    publication_date = DateTimeField()
    submission_id = StringField()
    last_update = DateTimeField()
    submission_date = DateTimeField()
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)


# Pubmed
class Pubmed(EmbeddedDocument):
    """The contents of a Pubmed document.

    This document contains specific information about publications.

    Attributes
    ----------
    pubmed_id: mongoengine.StringField
        The primary identifier for Pubmed. These are the accession number
        which begin with PMID.

    title: mongoengine.StringField
        Title of the paper.

    abstract: mongoengine.StringField
        Paper abstract.

    authors: mongoengine.ListField
        List of authors.

    citation: mongoengine.StringField
        Citation information for the paper.

    date_created: mongoengine.DateTimeField
        Date the pubmed entry was created.

    date_completed: mongoengine.DateTimeField
        Date the pubmed entry was completed.

    date_revised: mongoengine.DateTimeField
        Date the pubmed entry was last updated.

    """
    pubmed_id = StringField()
    title = StringField()
    abstract = StringField()
    authors = ListField(DictField())
    citation = StringField()
    date_created = DateTimeField()
    date_completed = DateTimeField()
    date_revised = DateTimeField()


# GEO
class Geo(EmbeddedDocument):
    """Subdocument for GEO."""
    pass


# Main Document class
class Ncbi(Document):
    """Document class that contains data from various NCBI databases."""
    srx = StringField(primary_key=True)
    sra = EmbeddedDocumentField(Sra)
    biosample = ListField(EmbeddedDocumentField(BioSample))
    bioproject = EmbeddedDocumentField(BioProject)
    geo = EmbeddedDocumentField(Geo)
    pubmed = ListField(EmbeddedDocumentField(Pubmed), default=list)

    meta = {'allow_inheritance': True}
