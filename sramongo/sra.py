"""Module to handle parsing SRA XML metadata."""

import re
from xml.etree import ElementTree
from collections import defaultdict
import datetime

from sramongo.logger import logger
from sramongo.xml_helpers import valid_path, parse_tree_from_dict
from sramongo.sra_const import EXISTING_STUDY_TYPES_ACTIVE, \
    EXISTING_STUDY_TYPES_DEPRICATED, INSTRUMENT_MODEL_ACTIVE, \
    INSTRUMENT_MODEL_DEPRICATED, LIBRARY_LAYOUT, LIBRARY_SELECTION, \
    LIBRARY_SOURCE, LIBRARY_STRATEGY, PLATFORMS


class XMLSchemaException(Exception):
    """Unexpected value given SRA xml schema definition."""
    pass


class SraExperiment(object):
    def __init__(self, node):
        """Class to parse SRA experiments and add them to mongoDB.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Experiment level node as an ElemenTree element.
        """
        # parse different sections for SRA xml
        self.sra = {
                'organization': self._parse_organization(node.find('Organization')),
                'submission': self._parse_submission(node.find('SUBMISSION')),
                'study': self._parse_study(node.find('STUDY')),
                'experiment': self._parse_experiment(node.find('EXPERIMENT')),
                'run': self._parse_run(node.find('RUN_SET')),
                'sample': self._parse_sample(node.find('SAMPLE')),
                'pool': self._parse_pool(node.find('Pool')),
                'db_flags': set(),
                'db_imported': datetime.datetime.now,
                }
        self.srx = self.sra['experiment']['experiment_id']

        # Add db flags
        if (self.sra['experiment']['library_source'] == 'TRANSCRIPTOMIC') | (self.sra['experiment']['library_strategy'] == 'RNA-Seq'):
            self.sra['db_flags'].add('RNASeq')

        if self.sra['experiment']['library_layout'] == 'SINGLE':
            self.sra['db_flags'].add('SE')
        elif self.sra['experiment']['library_layout'] == 'PAIRED':
            self.sra['db_flags'].add('PE')

    def _raise_xref_status(self, xref):
        """Some xrefs are more imporant and I want to pull them out.

        xrefs to specific databases I am interested in should be pulled out and
        stored separately. Go ahead and normalize database names and values a
        bit.

        Parameters
        ----------
        xref: dict
            A dictionary with keys 'db' and 'id'.
        """
        # databases of interest
        db_map = {'geo': 'GEO',
                  'gds': 'GEO_Dataset',
                  'bioproject': 'BioProject',
                  'biosample': 'BioSample',
                  'pubmed': 'pubmed'}

        # Make sure fully formed xref
        try:
            _id = xref['id']
            _db = xref['db']
        except:
            return False

        if (_id is None) | (_db is None):
            return False

        # normalize db name
        try:
            norm = _db.strip(' ()[].:').lower()
        except:
            norm = ''

        if norm in db_map.keys():
            # Normalize the ids a little
            id_norm = re.sub('geo|gds|bioproject|biosample|pubmed|pmid',
                             '', _id.lower()).strip(' :().').upper()
            return db_map[norm], id_norm
        else:
            return False

    def _drop_empty(self, d):
        """Scans through a dictionary and removes empy lists or None elements."""

        drop_keys = []
        for k, v in d.items():
            try:
                if (v is None) or (len(v) == 0):
                    drop_keys.append(k)
            except TypeError as err:
                # Some types (int, float) will cause an exception, if they are
                # of this type I don't really care then they are not blank so
                # skip.
                pass

        for k in drop_keys:
            del d[k]

        return d

    @valid_path
    def _parse_link_parts(self, node):
        """Parse the different parts of a link.

        In SRA a link can have a {label, db, id, query, url}. This function
        iterates over the different attributes and pulls out the bits of
        information.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is '{URL,XREF,ENTRE,DDBJ,ENA}_LINK'.
        """
        d = dict()
        for child in node:
            d[child.tag.lower()] = child.text
        return d

    @valid_path
    def _update_link(self, node, _type, d):
        """Given a type of link, search for it and return a dict."""

        def search(_type):
            """Returns ElemenTree search string."""
            return '*/' + _type.strip('s').upper()

        for l in node.findall(search(_type)):
            link = self._parse_link_parts(l)
            _raise = self._raise_xref_status(link)
            if _raise:
                d[_raise[0]] = _raise[1]
            else:
                d[_type].append(link)
        return d

    @valid_path
    def _parse_links(self, node):
        """Parse various types of links in the SRA.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is '*_LINKS'.
        """
        d = defaultdict(list)
        d.update(self._update_link(node, 'url_links', d))
        d.update(self._update_link(node, 'xref_links', d))
        d.update(self._update_link(node, 'entrez_links', d))
        d.update(self._update_link(node, 'ddbj_links', d))
        d.update(self._update_link(node, 'ena_links', d))
        return d

    @valid_path
    def _parse_ids(self, node, namespace):
        """Helper function to parse IDENTIFIERS sections.

        Automatically tries to pull out important external and secondary
        identifiers. Otherwises stores other identifiers in a list.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is 'IDENTIFIERS'.
        namespace: str
            Name to use for the primary ID.
        """
        d = defaultdict(list)

        for _id in node:
            if _id.tag == 'PRIMARY_ID':
                d[namespace + '_id'] = _id.text
            elif _id.tag == 'SUBMITTER_ID':
                xref = {'db': _id.get('namespace'), 'id': _id.text}
                d[_id.tag.lower()].append(xref)
            elif _id.tag == 'UUID':
                xref = _id.text
                d[_id.tag.lower()].append(xref)
            else:
                # Other types of ids (external, secondary).
                xref = {'db': _id.get('namespace'), 'id': _id.text}
                _raise = self._raise_xref_status(xref)
                if _raise:
                    d[_raise[0]] = _raise[1]
                else:
                    d[_id.tag.lower()].append(xref)
        return d

    # Organization
    @valid_path
    def _parse_organization(self, node):
        """Parses organization section."""
        d = dict()
        locs = {
            'organization_type': ('.', 'type'),
            'abbreviation': ('Name', 'abbr'),
            'name': ('Name', 'text'),
            'email': ('Contact', 'email'),
            'first_name': ('Contact/Name/First', 'text'),
            'last_name': ('Contact/Name/Last', 'text'),
        }
        d.update(parse_tree_from_dict(node, locs))
        d = self._drop_empty(d)
        return d

    # Submission
    @valid_path
    def _parse_submission(self, node):
        """Parses submission section."""
        d = dict()
        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'submission'))
        d.update(parse_tree_from_dict(node, {'broker': ('.', 'broker_name')}))
        d = self._drop_empty(d)
        return d

    @valid_path
    def _parse_study_type(self, node):
        """Processes study types section.

        Do XML related validation here based on the study type information and
        sra_const.

        node: xml.etree.ElementTree.ElementTree.element
            'STUDY/DESCRIPTOR/STUDY_TYPE'
        """
        d = dict()
        if node.get('existing_study_type'):
            d['study_type'] = node.get('existing_study_type')

            # XSD Validation
            # If depricated replace with active type
            d['study_type'] = EXISTING_STUDY_TYPES_DEPRICATED.get(d['study_type'], d['study_type'])

            # Make sure study type is in current active types
            if not d['study_type'] in EXISTING_STUDY_TYPES_ACTIVE:
                raise XMLSchemaException('Study type: {}'.format(d['study_type']))

        elif node.get('new_study_type'):
            d['study_type'] = node.get('new_study_type')

        return d

    @valid_path
    def _parse_study_related_study(self, node):
        """Parses related study section."""
        links = []
        for study in node:
            d = dict()
            d['db'] = study.find('RELATED_LINK/DB').text
            d['id'] = study.find('RELATED_LINK/ID').text
            d['label'] = study.find('RELATED_LINK/LABEL').text
            d['is_primary'] = study.find('IS_PRIMARY').text
            links.append(d)

        d = {'related_studies': links}
        d = self._drop_empty(d)
        return d

    # Study
    @valid_path
    def _parse_study(self, node):
        """Parses study section."""
        d = dict()
        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'study'))
        d.update(self._parse_links(node.find('STUDY_LINKS')))
        d.update(self._parse_study_type(node.find('DESCRIPTOR/STUDY_TYPE')))
        d.update(self._parse_study_related_study(
            node.find('DESCRIPTOR/RELATED_STUDIES')))

        locs = {
            'title': ('DESCRIPTOR/STUDY_TITLE', 'text'),
            'abstract': ('DESCRIPTOR/STUDY_ABSTRACT', 'text'),
            'center_name': ('.', 'center_name'),
            'center_project_name': (
                'DESCRIPTOR/CENTER_PROJECT_NAME', 'text'),
            'description': ('DESCRIPTOR/STUDY_DESCRIPTION', 'text'),
        }
        d.update(parse_tree_from_dict(node, locs))
        d = self._drop_empty(d)
        return d

    @valid_path(rettype=list)
    def _parse_attributes(self, node):
        """Parse various attributes as tag:value pairs.

        Parameters
        ----------
        node: xml.etree.ElementTree.ElementTree.Element
            Node where the root is 'EXPERIMENT_ATTRIBUTES'.
        """
        attributes = []
        for attribute in node.getchildren():
            attributes.append({'name': attribute.find('TAG').text,
                      'value': attribute.find('VALUE').text})
        return attributes

    # Sample
    @valid_path
    def _parse_sample(self, node):
        """Parses sample section."""
        d = dict()

        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'sample'))
        d.update(self._parse_links(node.find('SAMPLE_LINKS')))
        d['attributes'] = self._parse_attributes(node.find('SAMPLE_ATTRIBUTES'))

        locs = {
            'title': ('TITLE', 'text'),
            'taxon_id': ('SAMPLE_NAME/TAXON_ID', 'text'),
            'scientific_name': ('SAMPLE_NAME/SCIENTIFIC_NAME', 'text'),
            'common_name': ('SAMPLE_NAME/COMMON_NAME', 'text'),
            'individual_name': ('SAMPLE_NAME/INDIVIDUAL_NAME', 'text'),
            'description': ('SAMPLE/DESCRIPTION', 'text'),
        }
        d.update(parse_tree_from_dict(node, locs))
        d = self._drop_empty(d)
        return d

    @valid_path(rettype=list)
    def _parse_pool(self, node):
        """Parses pool section."""
        pool = []
        for member in node:
            d = dict()
            d.update(self._parse_ids(member.find('IDENTIFIERS'), 'sample'))
            pool.append(d['sample_id'])
        return pool

    # Experiment
    @valid_path
    def _parse_experiment(self, node):
        """Parses experiment section."""
        d = dict()

        d.update(self._parse_ids(node.find('IDENTIFIERS'), 'experiment'))
        d.update(self._parse_links(node.find('EXPERIMENT_LINKS')))
        d['attributes'] = self._parse_attributes(node.find('EXPERIMENT_ATTRIBUTES'))

        locs = {
            'title': ('TITLE', 'text'),
            'study_id': ('STUDY_REF/IDENTIFIERS/PRIMARY_ID', 'text'),
            'design': ('DESIGN/DESIGN_DESCRIPTION', 'text'),
            'library_name': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME', 'text'),
            'library_strategy': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY', 'text'),
            'library_source': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE', 'text'),
            'library_selection': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION', 'text'),
            'library_layout': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT',
                'child', 'tag'),
            'library_layout_orientation': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED',
                'ORIENTATION'),
            'library_layout_length': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED',
                'NOMINAL_LENGTH'),
            'library_layout_sdev': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED',
                'NOMINAL_SDEV'),
            'pooling_stategy': (
                'DESIGN/LIBRARY_DESCRIPTOR/POOLING_STRATEGY',
                'text'),
            'library_construction_protocol': (
                'DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL',
                'text'),
            'platform': ('PLATFORM', 'child', 'tag'),
            'instrument_model': ('PLATFORM/*/INSTRUMENT_MODEL', 'text'),
        }

        d.update(parse_tree_from_dict(node, locs))

        # XSD Validation
        def _clean_const(dat, name):
            consts = {'library_strategy': LIBRARY_STRATEGY,
                      'library_source': LIBRARY_SOURCE,
                      'library_selection': LIBRARY_SELECTION,
                      'library_layout': LIBRARY_LAYOUT,
                      'platform': PLATFORMS,
                      'instrument_model': INSTRUMENT_MODEL_ACTIVE,
                    }

            if dat[name].lower() in map(str.lower, consts[name]):
                for c in consts[name]:
                    if dat[name].lower() == c.lower():
                        return {name: c}
            else:
                logger.error('"{}" not found in {}.'.format(dat[name], name))
                return {}

        d.update(_clean_const(d, 'library_strategy'))
        d.update(_clean_const(d, 'library_source'))
        d.update(_clean_const(d, 'library_selection'))
        d.update(_clean_const(d, 'library_layout'))
        d.update(_clean_const(d, 'platform'))

        # Update instrument model if depricated
        d['instrument_model'] = INSTRUMENT_MODEL_DEPRICATED.get(d['instrument_model'], d['instrument_model'])
        d.update(_clean_const(d, 'instrument_model'))
        d = self._drop_empty(d)
        return d

    @valid_path
    def _parse_taxon(self, node):
        """Parse taxonomy informaiton."""

        def crawl(node, d=defaultdict(list)):
            for i in node:
                rank = i.get('rank', 'Unkown')
                d[rank].append({
                    'name': i.get('name').replace('.', '_').replace('$', ''),
                    'parent': node.get('name'),
                    'total_count': int(i.get('total_count')),
                    'self_count': int(i.get('self_count')),
                    'tax_id': i.get('tax_id'),
                    })
                if i.getchildren():
                    d.update(crawl(i, d))
            return d

        d = {'tax_analysis': {'nspot_analyze': node.get('analyzed_spot_count'),
                              'total_spots': node.get('total_spot_count'),
                              'mapped_spots': node.get('identified_spot_count'),
                              'tax_counts': crawl(node),
                              },
             }

        try:
            if d['tax_analysis']['nspot_analyze'] is not None:
                d['tax_analysis']['nspot_analyze'] = int(d['tax_analysis']['nspot_analyze'])
        except:
            logger.debug('Non integer count: nspot_analyze')
            logger.debug(d['tax_analysis']['nspot_analyze'])
            d['tax_analysis']['nspot_analyze'] = None

        try:
            if d['tax_analysis']['total_spots'] is not None:
                d['tax_analysis']['total_spots'] = int(d['tax_analysis']['total_spots'])
        except:
            logger.debug('Non integer count: total_spots')
            logger.debug(d['tax_analysis']['total_spots'])
            d['tax_analysis']['total_spots'] = None

        try:
            if d['tax_analysis']['mapped_spots'] is not None:
                d['tax_analysis']['mapped_spots'] = int(d['tax_analysis']['mapped_spots'])
        except:
            logger.debug('Non integer count: mapped_spots')
            logger.debug(d['tax_analysis']['mapped_spots'])
            d['tax_analysis']['mapped_spots'] = None

        return d

    @valid_path
    def _parse_run_reads(self, node, run_id):
        """Parse reads from runs."""
        d = dict()
        d['run_flags'] = set()

        nreads = node.get('nreads')
        try:
            if nreads is not None:
                d['nreads'] = int(nreads)
        except:
            # Sometimes nreads is set to 'variable' and no read information is
            # provided. I set these values to -1.
            if nreads == 'variable':
                d['nreads'] = -1
                d['run_flags'].add('no_read_information')
            else:
                raise ValueError('nreads is "{}", this value is not expected. '
                                 'Please adjust "sra.py" to account for this.'.format(nreads))

        if d['nreads'] == 1:
            # Single-ended Reads
            d['run_flags'].add('SE')
            read = node.find('Read')
            d['read_count_r1'] = float(read.get('count'))
            d['read_len_r1'] = float(read.get('average'))

        elif d['nreads'] == 2:
            # Pair-ended Reads
            d['read_count_r1'] = 0.0
            d['read_count_r2'] = 0.0
            d['read_len_r1'] = 0.0
            d['read_len_r2'] = 0.0

            rid = {"0": '_r1', "1": '_r2'}
            for read in node.findall('Read'):
                index = rid[read.get('index')]
                d['read_count' + index] = float(read.get('count'))
                d['read_len' + index] = float(read.get('average'))

            # Sometime nreads will be 2, but read len/count will be 0 for one
            # of the reads. Check for this and make SE if they are 0's
            if (d['read_len_r1'] != 0) & (d['read_len_r2'] != 0) & (d['read_count_r1'] != 0) & (d['read_count_r2'] != 0):
                d['run_flags'].add('PE')

                # Flag if read count or len are different
                if abs(d['read_len_r1'] - d['read_len_r2']) > 5:
                    d['run_flags'].add('pe_reads_not_equal_len')

                if abs(d['read_count_r1'] - d['read_count_r2']) > 10:
                    d['run_flags'].add('pe_reads_not_equal_count')
            else:
                d['run_flags'].add('SE')

                # Set which ever read has information to R1 and delete R2.
                if (d['read_len_r1'] == 0) | (d['read_count_r1'] == 0):
                    d['read_len_r1'] = d['read_len_r2']
                    d['read_count_r1'] = d['read_count_r2']
                    del d['read_len_r2']
                    del d['read_count_r2']

                elif (d['read_len_r2'] == 0) | (d['read_count_r2'] == 0):
                    del d['read_len_r2']
                    del d['read_count_r2']
        return d

    # Run
    @valid_path(rettype=list)
    def _parse_run(self, node):
        """Parses run section."""
        runs = []
        for run in node.findall('RUN'):
            d = dict()
            d.update(self._parse_ids(run.find('IDENTIFIERS'), 'run'))
            d['samples'] = self._parse_pool(run.find('Pool'))
            d.update(self._parse_taxon(run.find('tax_analysis')))
            d.update(self._parse_run_reads(run.find('Statistics'), d['run_id']))

            locs = {
                'experiment_id': ('EXPERIMENT_REF', 'accession'),
                'nspots': ('Statistics', 'nspots'),
                'nbases': ('Bases', 'count'),
            }

            d.update(parse_tree_from_dict(run, locs))

            # Make counts numeric
            if 'nspots' in d:
                try:
                    if d['nspots'] is not None:
                        d['nspots'] = int(d['nspots'])
                except:
                    logger.debug('Non numeric nspots: %s' % d['nspots'])
                    logger.debug('Setting nspots to None')
                    d['nspots'] = None

            if 'nbases' in d:
                try:
                    if d['nbases'] is not None:
                        d['nbases'] = int(d['nbases'])
                except:
                    logger.debug('Non numeric nbases: %s' % d['nbases'])
                    logger.debug('Setting nbases to None')
                    d['nbases'] = None

            d = self._drop_empty(d)
            runs.append(d)

        return runs
