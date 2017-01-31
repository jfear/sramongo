"""Set up MonoDB schema using monogodngine."""
from textwrap import fill
import datetime

from mongoengine import Document, EmbeddedDocument
from mongoengine import StringField, IntField, FloatField, \
    ListField, DictField, MapField, DateTimeField
from mongoengine import EmbeddedDocumentField, ReferenceField, signals
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

    def __str__(self):
        return DocumentString(self).string


class Xref(EmbeddedDocument):
    db = StringField()
    id = StringField()

    def __str__(self):
        return DocumentString(self).string


class XrefLink(EmbeddedDocument):
    label = StringField()
    db = StringField()
    id = StringField()
    meta = {'allow_inheritance': True}

    def __str__(self):
        return DocumentString(self).string


class EntrezLink(XrefLink):
    query = StringField()

    def __str__(self):
        return DocumentString(self).string


class DDBJLink(XrefLink):
    url = StringField()

    def __str__(self):
        return DocumentString(self).string


class ENALink(XrefLink):
    url = StringField()

    def __str__(self):
        return DocumentString(self).string


# Study
class Submission(EmbeddedDocument):
    """Submission embedded document.

    A submission must have a submission id (SRA/ERA/DRA). Additional metadata
    may be present about the submitter and external links to other databases.
    """
    submission_id = StringField(primary_key=True)
    broker = StringField()

    # Other IDS
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)
    uuid = ListField(StringField(), default=list)

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_SraExperiment(cls, sraExperiment, **kwargs):
        """Builds submission from an sramongo.SraExperiment.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no submission_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        sraExperiment: sramongo.SraExperiment
            An sra object parsed from XML.
        kwargs:
            Other name arguments will be used to update the sraExperiment prior
            to building.
        """
        try:
            sra = sraExperiment.submission
            sra.update(kwargs)
            submission = cls(**sra)
            return submission
        except ValidationError as err:
            logger.warn('%s\nSkipping this submission.' % err)
            logger.debug(sra)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sra)
            raise err


class Organization(EmbeddedDocument):
    """Organization embedded document.

    An organization has not defined id, so there are no required fields. This
    may result in multiple copies in the database.
    """
    organization_type = StringField()
    abbreviation = StringField()
    name = StringField()
    email = StringField()
    first_name = StringField()
    last_name = StringField()
    # TODO: Look into creating a primary key using first_name and last_name

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_SraExperiment(cls, sraExperiment, **kwargs):
        """Builds organization from an sramongo.SraExperiment.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        sraExperiment: sramongo.SraExperiment
            An sra object parsed from XML.
        kwargs:
            Other name arguments will be used to update the sraExperiment prior
            to building.
        """
        try:
            sra = sraExperiment.organization
            sra.update(kwargs)
            organization = cls(**sra)
            return organization
        except ValidationError as err:
            logger.warn('%s\nSkipping this organization.' % err)
            logger.debug(sra)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sra)
            raise err


class RelatedStudy(XrefLink):
    is_primary = StringField()


class Study(Document):
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

    submission: mongoengine.EmbeddedDocumentField
        A dictionary describing attributes about the Submission.

    organization = mongoengine.EmbeddedDocumentField
        A dictionary describing attributes about the Organization.

    experiments: mongoengine.ListField
        List of experiment_ids that are in this study.
    """
    study_id = StringField(primary_key=True)

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

    # NOTE: Additional Fields added post creation
    submission = EmbeddedDocumentField(Submission)
    organization = EmbeddedDocumentField(Organization)
    experiments = ListField(StringField(), default=list)

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_SraExperiment(cls, sraExperiment, **kwargs):
        """Builds study from an sramongo.SraExperiment.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no study_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        sraExperiment: sramongo.SraExperiment
            An sra object parsed from XML.
        kwargs:
            Other name arguments will be used to update the sraExperiment prior
            to building.
        """
        try:
            sra = sraExperiment.study
            sra.update(kwargs)
            if 'study_id' in sra:
                return cls.objects(pk=sra['study_id']).modify(upsert=True, new=True, **sra)
            else:
                raise ValidationError('No study_id')
        except ValidationError as err:
            logger.warn('%s\nSkipping this study.' % err)
            logger.debug(sra)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sra)
            raise err
        except Exception as err:
            logger.error(err)
            logger.debug(sra)
            raise err


# Samples
class Sample(Document):
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
    sample_id = StringField(primary_key=True)

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    BioSample = ReferenceField('BioSample')
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
    attributes = DictField()

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)

    meta = {'allow_inheritance': True}

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_SraExperiment(cls, sraExperiment, **kwargs):
        """Builds sample from an sramongo.SraExperiment.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no sample_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        sraExperiment: sramongo.SraExperiment
            An sra object parsed from XML.
        kwargs:
            Other name arguments will be used to update the sraExperiment prior
            to building.
        """
        try:
            sra = sraExperiment.sample
            sra.update(kwargs)

            if 'BioSample' in sra:
                sra['BioSample'] = BioSample.objects(pk=sra['BioSample']).modify(
                    upsert=True, new=True, biosample_id=sra['BioSample'])

            if 'sample_id' in sra:
                return cls.objects(pk=sra['sample_id']).modify(upsert=True, new=True, **sra)
            else:
                raise ValidationError('No sample_id')
        except ValidationError as err:
            logger.warn('%s\nSkipping this sample.' % err)
            logger.debug(sra)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sra)
            raise err


# Experiment
class Experiment(Document):
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

    attributes: mongoengine.DictField
        A set of key:value pairs describing the experiment.

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

    study: mongoengine.ReferenceField
        References the corresponding Study.

    samples: mongoengine.ListField
        List of references to corresponding samples.

    runs: mongoengine.ListField
        List of run_ids

    db_flags: mongoengine.ListField
        List of :ref:`database_flags`.

    """
    # SRX/DRX/ERX
    experiment_id = StringField(primary_key=True)

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
    attributes = DictField()

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)

    # NOTE: Additional Fields added post creation
    study = ReferenceField(Study)
    samples = ListField(ReferenceField(Sample), default=list)
    runs = ListField(StringField(), default=list)
    db_flags = ListField(StringField(), default=list)
    pipeline_flags = ListField(StringField(), default=list)

    meta = {'allow_inheritance': True}

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_SraExperiment(cls, sraExperiment, **kwargs):
        """Builds experiment from an sramongo.SraExperiment.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no experiment_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        sraExperiment: sramongo.SraExperiment
            An sra object parsed from XML.
        kwargs:
            Other name arguments will be used to update the sraExperiment prior
            to building.
        """
        try:
            sra = sraExperiment.experiment
            sra.update(kwargs)

            if ('study' in sra) and (isinstance(sra['study'], str)):
                study = Study.objects(pk=sra['study']).modify(
                            upsert=True, new=True, pk=sra['study'])

            if 'samples' in sra:
                ss = []
                for _s in sra['samples']:
                    if isinstance(_s, str):
                        ss.append(Sample.objects(pk=_s).modify(
                            upsert=True, new=True, pk=_s))

                sra['samples'] = ss

            if 'experiment_id' in sra:
                return cls.objects(pk=sra['experiment_id']).modify(
                        upsert=True, new=True, **sra)
            else:
                raise ValidationError('No experiment_id')
        except ValidationError as err:
            logger.warn('%s\nSkipping this experiment.' % err)
            logger.debug(sra)
            raise err
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sra)
            raise err


# Run
class TaxRecord(EmbeddedDocument):
    parent = StringField()
    total_count = IntField()
    self_count = IntField()
    tax_id = StringField()
    rank = StringField()

    def __str__(self):
        return DocumentString(self).string


class TaxAnalysis(EmbeddedDocument):
    nspot_analyze = IntField()
    total_spots = IntField()
    mapped_spots = IntField()
    tax_counts = MapField(EmbeddedDocumentField(TaxRecord))

    def __str__(self):
        return DocumentString(self).string


@update_modified.apply
class Run(Document):
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
        but it is stored here as a flat dictionary for easier access to
        different classes.

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

    experiment: mongoengine.ReferenceField
        This is an embedded reference to the corresponding experiment.

    release_date: mongoengine.StringField
        Release date of the Run.

    load_date: mongoengine.StringField
        Date the Run was uploaded.

    size_MB: mongoengine.IntField
        Size of the Run file.

    download_path: mongoengine.StringField
        Download path of the Run file.

    db_flags: mongoengine.ListField
        These are custom flags that I add. I have created these flags by
        interpreting specific information from data provided by the SRA. A list
        of :ref:`database_flags` and their descriptions.

    db_created: mongoengine.DateTimeField
        The date this document was added to the mongo database. In other words
        when you downloaded the data.

    db_modified: mongoengine.DateTimeField
        This is the date the document was last modified in the mongo database.

    """
    # SRR/DRR/ERR
    run_id = StringField(primary_key=True)

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

    # NOTE: Additional Fields added post creation
    experiment = ReferenceField(Experiment)
    release_date = StringField()
    load_date = StringField()
    size_MB = IntField()
    download_path = StringField()
    db_flags = ListField(StringField(), default=list)
    pipeline_flags = ListField(StringField(), default=list)
    db_created = DateTimeField(default=datetime.datetime.now)
    db_modified = DateTimeField()

    meta = {'allow_inheritance': True}

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_SraExperiment(cls, sraExperiment, runinfo, **kwargs):
        """Builds run from an sramongo.SraExperiment and a runinfo table.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no run_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        sraExperiment: sramongo.SraExperiment
            An sra object parsed from XML.
        runinfo: pandas.DataFrame
            A runinfo table imported as a data frame.
        kwargs:
            Other name arguments will be used to update the sraExperiment prior
            to building.
        """
        runs = []
        for sra in sraExperiment.run:
            try:
                if 'experiment' in sra:
                    sra['experiment'] = Experiment.objects(pk=sra['experiment']).modify(
                        upsert=True, new=True, pk=sra['experiment'])

                if 'run_id' in sra:
                    run = cls.objects(pk=sra['run_id']).modify(upsert=True, new=True, **sra)
                else:
                    raise ValidationError('No run_id')

                try:
                    if runinfo.notnull().loc[run.run_id, 'ReleaseDate']:
                        run.modify(release_date=runinfo.loc[run.run_id, 'ReleaseDate'])
                except KeyError:
                    pass

                try:
                    if runinfo.notnull().loc[run.run_id, 'LoadDate']:
                        run.modify(load_date=runinfo.loc[run.run_id, 'LoadDate'])
                except KeyError:
                    pass

                try:
                    if runinfo.notnull().loc[run.run_id, 'size_MB']:
                        run.modify(size_MB=int(runinfo.loc[run.run_id, 'size_MB']))
                except KeyError:
                    pass

                try:
                    if runinfo.notnull().loc[run.run_id, 'download_path']:
                        run.modify(download_path=runinfo.loc[run.run_id, 'download_path'])
                except KeyError:
                    pass

                try:
                    run.modfiy(experiment=Experiment.objects(experiment_id=run.experiment_id).first())
                except:
                    pass

                runs.append(run)

            except ValidationError as err:
                logger.warn('%s\nSkipping this run.' % err)
                logger.debug(sra)
                logger.debug(runinfo.loc[run.run_id, :])
            except FieldDoesNotExist as err:
                logger.error(err)
                logger.debug(sra)
                raise err

        return runs


class Contacts(EmbeddedDocument):
    email = StringField()
    first_name = StringField()
    last_name = StringField()

    def __str__(self):
        return DocumentString(self).string


# BioSample
class BioSample(Document):
    """The contents of a BioSample.

    BioSample is another database housed at NCBI which records sample metadata.
    This information should already be present in the SRA.Sample information,
    but to be safe we can pull into the BioSample for additional metadata.

    Attributes
    ----------

    biosample_id: mongoengine.StringField
        The primary identifier for a BioSample. These are the accession number
        which begin with SAM.

    biosample_secondary: mongoengine.StringField
        A secondary identifier for a BioSample. Identifiers begin with
        SAM.

    db_id: mongoengine.StringField
        This is BioSample's unique ID number.

    sample_id: mongoengine.StringField
        Unique id of the sample provided by the submitter.

    GEO: mongoengine.StringField
        GEO sample ID (GSM).

    title: mongoengine.StringField
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

    attributes: mongoengine.DictField
        A set of key:value pairs describing the sample. For example tissue:ovary
        or sex:female.

    """
    biosample_id = StringField(primary_key=True)
    biosample_secondary = StringField()
    db_id = StringField()
    sample_id = StringField()
    GEO = StringField()
    title = StringField()
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
    attributes = DictField()

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_BioSample(cls, bioSample, **kwargs):
        """Builds BioSample from an sramongo.biosample.BioSample.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no biosample_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        bioSample: sramongo.biosample.BioSample
            An biosample object parsed from XML.
        kwargs:
            Other name arguments will be used to update the BioSample prior
            to building.
        """
        try:
            bio = bioSample.biosample

            if 'biosample_id' in bio:
                bs = cls.objects(pk=bio['biosample_id']).modify(
                            upsert=True, new=True, **bio)
                # Make sure the samples are pointing at the right BioSample
                # object
                Sample.objects(BioSample=bs.db_id).update(BioSample=bs)

                # Delete BioSample document with wrong id
                BioSample.objects(pk=bs.db_id).delete()

                return bs
            else:
                raise ValidationError('No biosample_id')

        except ValidationError as err:
            logger.warn('%s\nSkipping this BioSample.' % err)
            logger.debug(bio)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(bio)
            raise err
        except Exception as err:
            logger.error(err)
            logger.debug(bio)
            raise err


# BioProject
class BioProject(Document):
    """The contents of a BioProject.

    BioProject is another database housed at NCBI which records project
    metadata.  This information should already be present in the SRA
    information, but to be safe we can pull into the BioProject for additional
    metadata.

    Attributes
    ----------
    bioproject_id: mongoengine.StringField
        The primary identifier for a BioProject. These are the accession number
        which begin with PRJ.

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
    bioproject_id = StringField(primary_key=True)
    name = StringField()
    title = StringField()
    description = StringField()
    publication = StringField()
    publication_date = DateTimeField()
    submission_id = StringField()
    last_date = DateTimeField()
    submission_date = DateTimeField()
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_BioProject(cls, bioProject, **kwargs):
        """Builds BioProject from an sramongo.bioproject.BioProject.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no bioproject_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        bioProject: sramongo.bioproject.BioProject
            An bioproject object parsed from XML.
        kwargs:
            Other name arguments will be used to update the BioProject prior
            to building.
        """
        try:
            bio = bioProject.bioproject
            if 'bioproject_id' in bio:
                return cls.objects(pk=bio['bioproject_id']).modify(
                            upsert=True, new=True, **bio)
            else:
                raise ValidationError('No bioproject_id')

        except ValidationError as err:
            logger.warn('%s\nSkipping this BioProject.' % err)
            logger.debug(bio)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(bio)
            raise err
        except Exception as err:
            logger.error(err)
            logger.debug(bio)
            raise err


# Pubmed
class Pubmed(Document):
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
    pubmed_id = StringField(primary_key=True)
    title = StringField()
    abstract = StringField()
    authors = ListField(DictField())
    citation = StringField()
    date_created = DateTimeField()
    date_completed = DateTimeField()
    date_revised = DateTimeField()

    def __str__(self):
        return DocumentString(self).string

    @classmethod
    def build_from_Pubmed(cls, pubmed, **kwargs):
        """Builds Pubmed from an sramongo.pubmed.Pubmed.

        Pulls in information and tries to validate. If there is a
        ValidationError (i.e. no bioproject_id or additional fields that have
        not been defined) then return None.

        Parameters
        ----------
        pubmed: sramongo.pubmed.Pubmed
            An pubmed object parsed from XML.
        kwargs:
            Other name arguments will be used to update the BioProject prior
            to building.
        """
        try:
            pub = Pubmed.pubmed
            if 'pubmed_id' in pub:
                return cls.objects(pk=pub['pubmed_id']).modify(
                            upsert=True, new=True, **pub)
            else:
                raise ValidationError('No pubmed_id')
        except ValidationError as err:
            logger.warn('%s\nSkipping this Pubmed.' % err)
            logger.debug(pub)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(pub)
            raise err
        except Exception as err:
            logger.error(err)
            logger.debug(bio)
            raise err
