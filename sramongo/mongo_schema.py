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
        document. In other words allows prettry printing of documents.

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
        representations. Carries indention level through other methods so that
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
            sraSubmission = sraExperiment.submission
            sraSubmission.update(kwargs)
            submission = cls(**sraSubmission)
            return submission
        except ValidationError as err:
            logger.warn('%s\nSkipping this submission.' % err)
            logger.debug(sraSubmission)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sraSubmission)
            raise err


class Organization(EmbeddedDocument):
    """Organization embedded document.

    An organization has not defined id, so there are no required fields. This
    may result in multiple copies in the database.
    """
    type = StringField()
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
            sraOrganization = sraExperiment.organization
            sraOrganization.update(kwargs)
            organization = cls(**sraOrganization)
            return organization
        except ValidationError as err:
            logger.warn('%s\nSkipping this organization.' % err)
            logger.debug(sraOrganization)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sraOrganization)
            raise err


class RelatedStudy(XrefLink):
    is_primary = StringField()


class Study(Document):
    """Study Document

    A study must have a study id (SRP/ERP/DRP). Additional metadata
    may be present. Studies embed submission and organization information.
    """
    # SRA/DRA/ERA
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
    type = StringField()
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
            sraStudy = sraExperiment.study
            sraStudy.update(kwargs)
            study = cls(**sraStudy)
            study.save()
            return study
        except ValidationError as err:
            logger.warn('%s\nSkipping this study.' % err)
            logger.debug(sraStudy)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sraStudy)
            raise err


# Samples
class Sample(Document):
    """Sample Document.

    A sample must have a sample id (SRS/ERS/DRS). Additional metadata
    may be present along with external links to other databases.
    """
    # SRS/DRS/ERS
    sample_id = StringField(primary_key=True)

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
    attributes = DictField()

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)

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
            sraSample = sraExperiment.sample
            sraSample.update(kwargs)
            sample = cls(**sraSample)
            sample.save()
            return sample
        except ValidationError as err:
            logger.warn('%s\nSkipping this sample.' % err)
            logger.debug(sraSample)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sraSample)
            raise err


# Experiment
class Experiment(Document):
    """Experiment Document.

    A experiment must have a experiment id (SRX/ERX/DRX). Additional metadata
    may be present along with external links to other databases.
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
            sraExp = sraExperiment.experiment
            sraExp.update(kwargs)
            experiment = cls(**sraExp)
            experiment.save()
            return experiment
        except ValidationError as err:
            logger.warn('%s\nSkipping this experiment.' % err)
            logger.debug(sraExp)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(sraExp)
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

    A run must have a run id (SRR/ERR/DRR). Additional metadata
    may be present about the submitter and external links to other databases.
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

    # if single ended
    read_count = FloatField()
    read_len = FloatField()

    # if paired ended
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
        for sraRun in sraExperiment.run:
            try:
                run = cls(**sraRun)

                try:
                    if runinfo.notnull().loc[run.run_id, 'ReleaseDate']:
                        run.release_date = runinfo.loc[run.run_id, 'ReleaseDate']
                except KeyError:
                    pass

                try:
                    if runinfo.notnull().loc[run.run_id, 'LoadDate']:
                        run.load_date = runinfo.loc[run.run_id, 'LoadDate']
                except KeyError:
                    pass

                try:
                    if runinfo.notnull().loc[run.run_id, 'size_MB']:
                        run.size_MB = int(runinfo.loc[run.run_id, 'size_MB'])
                except KeyError:
                    pass

                try:
                    if runinfo.notnull().loc[run.run_id, 'download_path']:
                        run.download_path = runinfo.loc[run.run_id, 'download_path']
                except KeyError:
                    pass

                try:
                    run.experiment = Experiment.objects(experiment_id=run.experiment_id).first()
                except:
                    pass

                run.save()
                runs.append(run)

            except ValidationError as err:
                logger.warn('%s\nSkipping this run.' % err)
                logger.debug(sraRun)
                logger.debug(runinfo.loc[run.run_id, :])
            except FieldDoesNotExist as err:
                logger.error(err)
                logger.debug(sraRun)
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
    """BioSample Document.

    A BioSample must have a run id (SAMN). Additional metadata may be present.
    """
    biosample_id = StringField(primary_key=True)
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
    contacts = ListField(EmbeddedDocumentField(Contacts))
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
            biosample = cls(**bioSample)
            biosample.save()
        except ValidationError as err:
            logger.warn('%s\nSkipping this BioSample.' % err)
            logger.debug(bioSample)
            return None
        except FieldDoesNotExist as err:
            logger.error(err)
            logger.debug(bioSample)
            raise err

        return biosample
