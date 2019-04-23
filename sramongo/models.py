"""Document models for MongoDB"""
from datetime import datetime

from mongoengine import Document, EmbeddedDocument
from mongoengine import StringField, IntField, FloatField, ListField, DictField, DateTimeField, MapField
from mongoengine import EmbeddedDocumentField


class Attribute(EmbeddedDocument):
    name = StringField()
    value = StringField()


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


class Study(EmbeddedDocument):
    """The contents of a SRA study.

    A study consists of a set of experiments designed with an overall goal in
    mind. For example, this could include a control experiment and a treatment
    experiment with the goal being to identify expression differences resulting
    from the treatment. The SRA study is the top level of the submission
    hierarchy.

    Attributes
    ----------
    accn: mongoengine.StringField
        The primary identifier for a study. Identifiers begin with
        SRP/ERP/DRP depending on which database they originate from.

    bioproject: mongoengine.StringField
        The associated BioProject identifier.

    geo: mongoengine.StringField
        The associated GEO identifier.

    geo: mongoengine.StringField
        The associated Pubmed identifiers.

    title: mongoengine.StringField
        The title of the study.

    abstract: mongoengine.StringField
        Abstract of the study.

    center_name: mongoengine.StringField
        Name of the submitting center.

    center_project_name: mongoengine.StringField
        Center specific identifier for the study.

    description: mongoengine.StringField
        Additional text describing the study.

    """
    accn = StringField()

    # External IDs
    bioproject = StringField()
    geo = StringField()
    pubmed = ListField(IntField())

    # Attributes
    title = StringField()
    abstract = StringField()
    center_name = StringField()
    center_project_name = StringField()
    description = StringField()


class Sample(EmbeddedDocument):
    """The contents of a SRA sample.

    A sample is the biological unit. An individual sample or a pool of samples
    can be use in the SRA Experiment. This document contains information
    describing the sample ranging from species information to detailed
    descriptions of what and how material was collected.

    Attributes
    ----------
    accn: mongoengine.StringField
        The primary identifier for a sample. Identifiers begin with
        SRS/ERS/DRS depending on which database they originate from.

    biosample: mongoengine.StringField
        The associated BioSample identifier.

    geo: mongoengine.StringField
        The associated GEO identifier.

    title: mongoengine.StringField
        The title of the sample.

    taxon_id: mongoengine.IntField
        The NCBI taxon id.

    scientific_name: mongoengine.StringField
        The scientific name.

    common_name: mongoengine.StringField
        The common name.

    attributes: mongoengine.DictField
        A set of key:value pairs describing the sample. For example tissue:ovary
        or sex:female.

    """
    # SRS/DRS/ERS
    accn = StringField()

    # External IDs
    biosample = StringField()
    geo = StringField()

    # Attributes
    title = StringField()
    taxon_id = IntField()
    scientific_name = StringField()
    common_name = StringField()
    attributes = ListField(EmbeddedDocumentField(Attribute), default=list)


class Run(EmbeddedDocument):
    """Run Document.

    A Run describes a dataset generated from an Experiment. For example if a
    Experiment is sequenced on multiple lanes of a Illumina flowcell then data
    from each lane are considered a Run.

    Attributes
    ----------
    srr: mongoengine.StringField
        The primary identifier for a run. Identifiers begin with
        SRR/ERR/DRR depending on which database they originate from.

    nspots: mongoengine.IntField
        The total number of spots on a Illumina flowcell.

    nbases: mongoengine.IntField
        The number of bases.

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

    release_date: mongoengine.DateTimeField
        Release date of the Run. This information is from the runinfo table and
        not the XML.

    load_date: mongoengine.DateTimeField
        Date the Run was uploaded. This information is from the runinfo table
        and not the XML.

    size_MB: mongoengine.IntField
        Size of the Run file. This information is from the runinfo table and not
        the XML.

    """
    # SRR/DRR/ERR
    srr = StringField()

    # Attributes
    nspots = IntField()
    nbases = IntField()
    nreads = IntField()

    # if single ended then just use _r1
    read_count_r1 = FloatField()
    read_len_r1 = FloatField()

    read_count_r2 = FloatField()
    read_len_r2 = FloatField()

    # NOTE: Additional Fields not in the SRA XML but in summary table
    release_date = DateTimeField()
    load_date = DateTimeField()
    size_MB = IntField()


class Geo(EmbeddedDocument):
    # TODO add geo parser
    accn = StringField()
    GEO_Dataset = StringField()
    sramongo_last_updated = DateTimeField(default=datetime.now())


class BioProject(EmbeddedDocument):
    """The contents of a BioProject.

    BioProject is another database housed at NCBI which records project
    metadata.  This information should already be present in the SRA
    information, but to be safe we can pull into the BioProject for additional
    metadata.

    Attributes
    ----------
    accn: mongoengine.StringField
        The primary identifier for a BioProject. These are the accession number
        which begin with PRJ.

    id: mongoengine.IntField
        The primary identifier for a BioProject. These are the id numbers.

    name: mongoengine.StringField
        A brief name of the project.

    title: mongoengine.StringField
        The title of the project.

    description: mongoengine.StringField
        A short description of the project.

    last_date: mongoengine.DateTimeField
        Last date the BioProject was updated.

    submission_date: mongoengine.DateTimeField
        Date the BioProject was submitted.

    """
    accn = StringField()
    bioproject_id = IntField()
    name = StringField()
    title = StringField()
    description = StringField()
    last_update = DateTimeField()
    submission_date = DateTimeField()

    sramongo_last_updated = DateTimeField(default=datetime.now())


class BioSample(EmbeddedDocument):
    """The contents of a BioSample.

    BioSample is another database housed at NCBI which records sample metadata.
    This information should already be present in the Sra.sample information,
    but to be safe we can pull into the BioSample for additional metadata.

    Attributes
    ----------
    accn: mongoengine.StringField
        The primary identifier for a BioSample. These are the accession number
        which begin with SAM.

    id: mongoengine.IntField
        The primary identifier for a BioSample. These are the id number.

    title: mongoengine.StringField
        A free text description of the sample.

    description: mongoengine.StringField
        A free text description of the sample.

    publication_date: mongoengine.StringField
        Date the sample was published.

    last_update: mongoengine.StringField
        Last time BioSample updated sample information.

    submission_date: mongoengine.StringField
        Date the sample was submitted

    attributes: mongoengine.ListField of mongoengine.DictField
        A list of dictionaries containing key:value pairs describing the
        experiment. The stored dictionaries are of the form {'name': value,
        'value': value}. This was done to make querying easier.

    """
    accn = StringField()
    biosample_id = IntField()
    title = StringField()
    description = StringField()
    last_update = StringField()
    submission_date = StringField()
    contacts = ListField(DictField(), default=list)
    attributes = ListField(EmbeddedDocumentField(Attribute), default=list)
    sramongo_last_updated = DateTimeField(default=datetime.now())


class Pubmed(EmbeddedDocument):
    """The contents of a Pubmed document.

    This document contains specific information about publications.

    Attributes
    ----------
    accn: mongoengine.StringField
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
    accn = StringField()
    title = StringField()
    abstract = StringField()
    authors = ListField(DictField())
    citation = StringField()
    date_created = DateTimeField()
    date_completed = DateTimeField()
    date_revised = DateTimeField()
    sramongo_last_updated = DateTimeField(default=datetime.now())


class SraDocument(Document):
    srx = StringField()
    sra_id = IntField()
    title = StringField()
    design = StringField()

    sramongo_last_updated = DateTimeField(default=datetime.now())

    sra_create_date = DateTimeField()
    sra_update_date = DateTimeField()

    # Technical Attributes
    library_name = StringField()
    library_strategy = StringField()
    library_source = StringField()
    library_selection = StringField()
    library_layout = StringField()
    library_layout_length = StringField()
    library_layout_sdev = StringField()
    library_construction_protocol = StringField()
    platform = StringField()
    instrument_model = StringField()

    # Embedded Documents
    organization = EmbeddedDocumentField(Organization)
    study = EmbeddedDocumentField(Study)
    sample = EmbeddedDocumentField(Sample)
    runs = ListField(EmbeddedDocumentField(Run))
    BioProject = EmbeddedDocumentField(BioProject)
    BioSmaple = EmbeddedDocumentField(BioSample)
    papers = ListField(EmbeddedDocumentField(Pubmed))
    Geo = EmbeddedDocumentField(Geo)


class TaxRecord(EmbeddedDocument):
    parent = StringField()
    total_count = IntField()
    self_count = IntField()
    tax_id = StringField()
    name = StringField()


class TaxAnalysis(Document):
    """
        A dictionary containing results from a taxonomic analysis. Some Runs are
        analyzed and the number of reads that align to different taxa are
        recorded. The taxanomic analysis is stored in the SRA as a hierarchy,
        but it is stored here as a flattend dictionary for easier access to
        different classes. Basic structure is:

        'nspoot_analyze': The number of spots analyzed,
        'total_spots': The total number of spots,
        'mapped_spots': The number of spots that were able to be mapped,
        'tax_count': A dictionary containing actual taxonomic counts organized by level in the tree of life

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
    """
    srr = StringField()
    nspot_analyze = IntField()
    total_spots = IntField()
    mapped_spots = IntField()
    tax_counts = MapField(ListField(EmbeddedDocumentField(TaxRecord), default=list))
