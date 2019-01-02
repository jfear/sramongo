"""Module for interacting with NCBI's ENTREZ API.

NCBI provides an API for querying and downloading data from their databases.

"""
import time
from typing import Optional, List
import urllib.parse
import requests
from collections import namedtuple
import re

from dateutil.parser import parse

BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

EsearchResult = namedtuple('EsearchResult', 'ids count webenv query_key')
EsummaryResult = namedtuple('EsummaryResult', 'id srx create_date update_date')


def esearch(database, query, userhistory=True, webenv=False, query_key=False, retstart=False, retmax=False) \
        -> Optional[EsearchResult]:
    """Search for a query using the Entrez ESearch API.

    Parameters
    ----------
    database : str
        Entez database to search.
    query : str
        Query string
    userhistory : bool
        Tells API to return a WebEnV and query_key.
    webenv : str
        An Entrez WebEnv to use saved history.
    query_key : str
        An Entrez query_key to use saved history.
    retstart : int
        Return values starting at this index.
    retmax : int
        Return at most this number of values.

    Returns
    -------
    EsearchResult
        A named tuple with values [ids, count, webenv, query_key]

    """
    cleaned_query = urllib.parse.quote_plus(query)

    url = BASE_URL + f'esearch.fcgi?db={database}&term={cleaned_query}&retmode=json'

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
        int(text['esearchresult'].get('count', '')),
        text['esearchresult'].get('webenv', ''),
        text['esearchresult'].get('querykey', '')
    )


def epost(database, ids: List[str], webenv=False) -> Optional[requests.Response]:
    """Post IDs using the Entrez ESearch API.

    Parameters
    ----------
    database : str
        Entez database to search.
    ids : list
        List of IDs to submit to the server.
    webenv : str
        An Entrez WebEnv to post ids to.

    Returns
    -------
    requests.Response

    """
    id = ','.join(ids)
    url_params = f'db={database}&id={id}&retmode=json'

    if webenv:
        url_params += f'&WebEnv={webenv}'

    if len(ids) <= 500:
        url = BASE_URL + f'epost.fcgi' + '?' + url_params
        resp = requests.get(url)
        if resp.status_code != 200:
            print('There was a server error')
            return
    else:
        url = BASE_URL + f'epost.fcgi'
        resp = requests.post(url, url_params)

    return resp


def esummary(database, ids=False, webenv=False, query_key=False, count=False, retstart=False, retmax=False) \
        -> Optional[List[EsummaryResult]]:
    """Get document summaries using the Entrez ESearch API.

    Parameters
    ----------
    database : str
        Entez database to search.
    ids : list or str
        List of IDs to submit to the server.
    webenv : str
        An Entrez WebEnv to use saved history.
    query_key : str
        An Entrez query_key to use saved history.
    count : int
        Number of records in the webenv
    retstart : int
        Return values starting at this index.
    retmax : int
        Return at most this number of values.

    Returns
    -------
    list
        A list of EsummaryResults with values [id, srx, create_date, update_date]

    """
    url = BASE_URL + f'esummary.fcgi?db={database}&retmode=json'

    if webenv and query_key:
        url += f'&WebEnv={webenv}&query_key={query_key}'
    elif ids:
        if isinstance(ids, str):
            id = ids
        else:
            id = ','.join(ids)
        url += f'&id={id}'
        count = len(id.split(','))

    results = []
    for resp in entrez_sets_of_results(url, retstart, retmax, count):
        text = resp.json()
        results.extend(parse_esummary(text))

    return results


def entrez_sets_of_results(url, retstart=False, retmax=False, count=False):
    """Gets sets of results back from Entrez.

    Entrez can only return 500 results at a time. This creates a generator that gets results by incrementing
    retstart and retmax.

    Parameters
    ----------
    url : str
        The Entrez API url to use.
    retstart : int
        Return values starting at this index.
    retmax : int
        Return at most this number of values.
    count : int
        The number of results returned by EQuery.

    Yields
    ------
    requests.Response

    """
    if not retstart:
        retstart = 0

    if not retmax:
        retmax = 500

    if not count:
        count = retmax

    retmax = 500  # Entrez can return a max of 500
    while retstart < count:
        diff = count - retstart
        if diff < 500:
            retmax = diff

        _url = url + f'&retstart={retstart}&retmax={retmax}'
        resp = esummary_try_three_times(_url)
        if resp is None:
            return

        retstart += retmax
        yield resp


def esummary_try_three_times(url, attempt=0, num_tries=3) -> Optional[requests.Response]:
    while attempt < num_tries:
        resp = requests.get(url)
        if resp.status_code != 200:
            attempt += 1
            time.sleep(1)
            esummary_try_three_times(url, attempt, num_tries)

        return resp

    print('There was a server error')
    return


def parse_esummary(json: dict) -> List[EsummaryResult]:
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


def efetch(database, ids=False, webenv=False, query_key=False, count=False, retstart=False, retmax=False,
           rettype='full', retmode='xml') -> Optional[List[EsummaryResult]]:
    """Get documents using the Entrez ESearch API.

    Parameters
    ----------
    database : str
        Entez database to search.
    ids : list or str
        List of IDs to submit to the server.
    webenv : str
        An Entrez WebEnv to use saved history.
    query_key : str
        An Entrez query_key to use saved history.
    count : int
        Number of records in the webenv
    retstart : int
        Return values starting at this index.
    retmax : int
        Return at most this number of values.
    rettype : str
        The type of document to return. Refer to link for valid return types for each database.
        https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
    retmode : str
        The format of document to return. Refer to link for valid formats for each database.
        https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

    Yields
    ------
    str
        Text from effect results. Format depends on parameters passed to retmode

    """
    url = BASE_URL + f'efetch.fcgi?db={database}&retmode={retmode}&rettype={rettype}'

    if webenv and query_key:
        url += f'&WebEnv={webenv}&query_key={query_key}'
    elif ids:
        if isinstance(ids, str):
            id = ids
        else:
            id = ','.join(ids)
        url += f'&id={id}'
        count = len(id.split(','))

    for resp in entrez_sets_of_results(url, retstart, retmax, count):
        yield resp.text

# TODO: Maybe add XML parsing here to make it look like a JSON object
