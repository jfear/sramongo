"""Set up MonoDB schema using monogodngine."""
from textwrap import fill
from mongoengine import Document, EmbeddedDocument
from mongoengine import StringField, ListField, DictField, MapField
from mongoengine import EmbeddedDocumentField, ReferenceField


class DocumentString(object):
    def __init__(self, document):
        """Represents mongo Document as strings.

        This is a helper class to make a string representation of a mongo
        document.

        Parameters
        ----------
        document: mongoengine.Document
            A subclass of the mongoengine.Document. Assumes that all data is
            stored in an attribute self._data.

        """
        self.string = ''
        self.add_dict(document._data)

    def add_dict(self, dictionary, spacer=0):
        """Iterates over a dictionary and parses strings.

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
        """Add key value pair to string for pretty printing."""
        spacer = ' ' * spacer
        if value is not None:
            # wrap really long lines like abstracts
            if len(value) > 100:
                s = spacer + ' ' * (len(key) + 2)
                value = fill(value, width=75).replace('\n', '\n' + s)
        self.string += spacer + '{}: {}\n'.format(key, value)

    def add_list(self, key, value, spacer=0):
        """Formats ids into astring for printing.

        Parameters
        ----------
        ids: list
            List of ids with a proper string representation (i.e. Xref)
        name: str
            Name of the id will be appended as a title (i.e. external_id)
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
        """Adds embed document."""
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
    submission_id = StringField(primary_key=True)
    broker = StringField()

    # Other IDS
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)

    def __str__(self):
        return DocumentString(self).string


class Organization(EmbeddedDocument):
    type = StringField()
    abbreviation = StringField()
    name = StringField()
    email = StringField()
    first_name = StringField()
    last_name = StringField()

    def __str__(self):
        return DocumentString(self).string


class RelatedStudy(XrefLink):
    is_primary = StringField()


class Study(Document):
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


# Samples
class Sample(Document):
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


# Experiment
class Pool(EmbeddedDocument):
    sample_id = ReferenceField(Sample)
    GEO = StringField()
    BioSample = StringField()

    def __str__(self):
        return DocumentString(self).string


class Experiment(Document):
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
    library_layout_nominal_length = StringField()
    library_layout_nominal_sdev = StringField()
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

    # NOTE: Additional Fields added post creation
    samples = ListField(EmbeddedDocumentField(Pool))

    # NOTE: Additional Fields added post creation
    runs = ListField(StringField(), default=list)

    def __str__(self):
        return DocumentString(self).string


# Run
class TaxRecord(EmbeddedDocument):
    parent = StringField()
    total_count = StringField()
    self_count = StringField()
    tax_id = StringField()
    rank = StringField()

    def __str__(self):
        return DocumentString(self).string


class TaxAnalysis(EmbeddedDocument):
    nspot_analyze = StringField()
    total_spots = StringField()
    mapped_spots = StringField()
    tax_counts = MapField(EmbeddedDocumentField(TaxRecord))

    def __str__(self):
        return DocumentString(self).string


class Run(Document):
    # SRR/DRR/ERR
    run_id = StringField(primary_key=True)

    # Other IDs
    external_id = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_id = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_id = ListField(EmbeddedDocumentField(Xref), default=list)

    # Attributes
    experiment_id = StringField()
    samples = ListField(EmbeddedDocumentField(Pool), default=list)
    nspots = StringField()
    nbases = StringField()
    tax_analysis = EmbeddedDocumentField(TaxAnalysis)
    nreads = StringField()

    # if single ended
    read_count = StringField()
    read_len = StringField()

    # if paired ended
    read_count_r1 = StringField()
    read_len_r1 = StringField()

    read_count_r2 = StringField()
    read_len_r2 = StringField()

    # NOTE: Additional Fields added post creation
    experiment = ReferenceField(Experiment)

    def __str__(self):
        return DocumentString(self).string
