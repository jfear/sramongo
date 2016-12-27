"""Set up MonoDB schema using monogodngine."""
from mongoengine import Document, EmbeddedDocument
from mongoengine import StringField, ListField, DictField
from mongoengine import EmbeddedDocumentField, ReferenceField


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


class TaxRecord(EmbeddedDocument):
    parent = StringField()
    total_count = StringField()
    self_count = StringField()
    tax_id = StringField()
    rank = StringField()


class TaxAnalysis(EmbeddedDocument):
    nspot_analyze = StringField()
    total_spots = StringField()
    mapped_spots = StringField()
    tax_counts = EmbeddedDocumentField(TaxRecord)


class RelatedStudy(XrefLink):
    is_primary = StringField()


class Submission(EmbeddedDocument):
    submission_id = StringField(required=True, unique=True)

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    GEO_Dataset = StringField()
    BioProject = StringField()
    BioSample = StringField()
    pubmed = StringField()

    # Other IDS
    external_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_ids = ListField(EmbeddedDocumentField(Xref), default=list)


class Organization(EmbeddedDocument):
    type = StringField()
    abbreviation = StringField()
    name = StringField()
    email = StringField()
    first_name = StringField()
    last_name = StringField()


class Study(Document):
    study_id = StringField(required=True, unique=True)

    # Looked through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    GEO_Dataset = StringField()
    BioProject = StringField()
    BioSample = StringField()
    pubmed = StringField()

    # Other IDs
    external_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_ids = ListField(EmbeddedDocumentField(Xref), default=list)

    # Attributes
    title = StringField()
    type = StringField()
    abstract = StringField()
    center_name = StringField()
    center_project_name = StringField()
    description = StringField()

    # Related study info
    related_studies = ListField(
            EmbeddedDocumentField(RelatedStudy), default=list)

    # External links
    url_links = ListField(EmbeddedDocumentField(URLLink), default=list)
    xref_links = ListField(EmbeddedDocumentField(XrefLink), default=list)
    entrez_links = ListField(EmbeddedDocumentField(EntrezLink), default=list)
    ddbj_links = ListField(EmbeddedDocumentField(DDBJLink), default=list)
    ena_links = ListField(EmbeddedDocumentField(ENALink), default=list)

    # NOTE: Additional Fields added post creation
    submission = EmbeddedDocumentField(Submission)

    # NOTE: Additional Fields added post creation
    organization = EmbeddedDocumentField(Organization)


class Pool(EmbeddedDocument):
    sample_id = StringField()
    GEO = StringField()
    BioSample = StringField()


class Experiment(Document):
    experiment_id = StringField(required=True, unique=True)

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    GEO_Dataset = StringField()
    BioProject = StringField()
    BioSample = StringField()
    pubmed = StringField()

    # Other IDs
    external_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_ids = ListField(EmbeddedDocumentField(Xref), default=list)

    # Attributes
    title = StringField()
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
    library_contruction_protocol = StringField('library_construction_protocol')
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


class Sample(Document):
    sample_id = StringField(required=True, unique=True)

    # Look through the external/secondary/submitter for these database xrefs
    # and pull them out for easy access.
    GEO = StringField()
    GEO_Dataset = StringField()
    BioProject = StringField()
    BioSample = StringField()
    pubmed = StringField()

    # Other IDs
    external_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_ids = ListField(EmbeddedDocumentField(Xref), default=list)

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


class Run(Document):
    run_id = StringField(required=True, unique=True)

    # Other IDs
    external_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    secondary_ids = ListField(EmbeddedDocumentField(Xref), default=list)
    submitter_ids = ListField(EmbeddedDocumentField(Xref), default=list)

    # Attributes
    experiment_id = ReferenceField(Experiment)
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


if __name__ == '__main__':
    pass
