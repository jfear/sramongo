"""Module for interacting with NCBI's ENTREZ API.

NCBI provides an API for querying and downloading data from their databases.

"""
import urllib.parse
import requests
from collections import namedtuple
import re
import typing

from dateutil.parser import parse

base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

EsearchResult = namedtuple('EsearchResult', 'ids webenv query_key')
EsummaryResult = namedtuple('EsummaryResult', 'id srx create_date update_date')

def esearch(database, query, userhistory=True, webenv=False, query_key=False, retstart=False, retmax=False):
    cleaned_query = urllib.parse.quote_plus(query)

    url = base_url + f'esearch.fcgi?db={database}&term={cleaned_query}&retmode=json'

    if userhistory:
        url += '&usehistory=y'

    if webenv:
        url += f'&WebEnv={webenv}'

    if query_key:
        url += f'&query_key={query_key}'

    if retstart:
        url += f'&retstart={retstart}'

    if retmax:
        url += f'&retmax={retmax}'

    resp = requests.get(url)
    if resp.status_code != 200:
        print('There was a server error')
        return

    text = resp.json()

    return EsearchResult(
        text['esearchresult'].get('idlist', []),
        text['esearchresult'].get('webenv', ''),
        text['esearchresult'].get('querykey', '')
    )


def epost(database, ids, webenv=False):
    id = ','.join(ids)
    url_params = f'db={database}&id={id}&retmode=json'

    if webenv:
        url_params += f'&WebEnv={webenv}'

    if len(ids) <= 500:
        url = base_url + f'epost.fcgi' + '?' + url_params
        resp = requests.get(url)
        if resp.status_code != 200:
            print('There was a server error')
            return
    else:
        url = base_url + f'epost.fcgi'
        resp = requests.post(url, url_params)

    return resp


def esummary(database, ids=False, webenv=False, query_key=False, retstart=False, retmax=False):

    # TODO: Need to add loop here because I can only get 500 documents back at a time.
    url = base_url + f'esummary.fcgi?db={database}&retmode=json'

    if webenv and query_key:
        url += f'&WebEnv={webenv}&query_key={query_key}'
    elif ids:
        if isinstance(ids, str):
            id = ids
        else:
            id = ','.join(ids)
        url += f'&id={id}'

    resp = requests.get(url)
    if resp.status_code != 200:
        print('There was a server error')
        return

    text = resp.json()
    return parse_esummary(text)


def parse_esummary(json) -> typing.List[EsummaryResult]:
    uids = json['result']['uids']

    srx_pattern = re.compile(r';Experiment acc=\"([SED]RX\d+)\"')

    results = []
    for uid in uids:
        xml: str = json['result'][uid].get('expxml', '')
        srx = re.findall(srx_pattern, xml)[0]
        create_date = parse(json['result'][uid].get('createdate', ''))
        update_date = parse(json['result'][uid].get('updatedate', ''))
        results.append(EsummaryResult(uid, srx, create_date, update_date))

    return results



def efetch():
    pass

