#!/usr/bin/env python3.6


# Copyright 2022 Robert Bossy (INRAE)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.


'''Wrapper objects and functions for the Omnicrobe REST API. See: http://omnicrobe.migale.inrae.fr/api.

Requirements:
    Python >= 3.6
    requests module, see: https://requests.readthedocs.io/

Variables:
    BASE_URL (str): default base URL of the Omnicrobe API.

'''

import requests


BASE_URL = 'https://omnicrobe.migale.inrae.fr/api'


class OmnicrobeError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


def fetch_json(fun_url, base_url=None, params=None):
    if base_url is None:
        base_url = BASE_URL
    url = base_url + fun_url
    resp = requests.get(url, params=params, verify=False)
    if resp.status_code == 200:
        return resp.json()
    raise OmnicrobeError(f'{resp.status_code} : {url}')


def version(base_url=None):
    '''Get the database version.

    Parameters:
        base_url (str): base URL of the API. Defaults to None. When this parameter is None, then use this module BASE_URL value.

    Returns:
        A dict object that contains the version of the API and each source in the database.
    '''
    return fetch_json('/get/version', base_url)


class LazyObject(object):
    '''A LazyObject holds an identifier and downloads other properties from the API when requested.

    Description:
        LazyObject is a class meant to be extended. Each Florolege entity type is represented as a class extending LazyObject.
        Override the methods _object_url() and _fill().

    Attributes:
        identifier (str): object identifier.
        base_url (str): base URL of the API.
        _loaded (boolean): either additional attributes have been loaded from the API.

        Other attributes may be filled from the API when accessed.
    '''

    def __init__(self, identifier, base_url=None):
        '''Creates a LazyObject object.

        Parameters:
            identifier (str): object identifier.
            base_url (str): base URL of the API. Defaults to None. When this parameter is None, then use this module BASE_URL value.
        '''
        self.identifier = identifier
        self.base_url = base_url
        self._loaded = False

    def __getattr__(self, name):
        '''Last resort object attribute access. If self._loaded is false, then this method attempts to load attributes from the API. See: https://docs.python.org/3/reference/datamodel.html#object.__getattr__ .'''
        if self._loaded:
            raise AttributeError(name)
        self._load()
        return getattr(self, name)

    def _load(self):
        '''
        Fetch from the Omnicrobe API thet data corresponding to this object. Calls self._fill() and set self._loaded to True.
        '''
        json = fetch_json(self._object_url(), self.base_url)
        self._fill(json)
        self._loaded = True

    def _object_url(self):
        '''Returns the URL in the Omnicrobe API corresponding to this object.'''
        raise NotImplementedError()

    def _fill(self, json):
        '''Modifies this object using the specified json object obtained from the Omnicrobe API.'''
        raise NotImplementedError()

    def convert_path(self, pathstr):
        '''Converts a slash-separated path of identifiers into a sequence of ancestor objects.

        Parameters:
            pathstr (str): identifier path sperated by slashes (/).

        Returns:
            A tuple of objects of the same type as this object.
        '''
        return tuple(self.__class__(identifier, self.base_url) for identifier in pathstr.split('/') if identifier != '')

    @staticmethod
    def path_string(path):
        return f'/{"/".join(a.identifier for a in path)}'


class Taxon(LazyObject):
    '''A Taxon object represents a taxon in the Omnicrobe database

    Attributes:
        Inherits from LazyObject.
        When the data is loaded from the API, then the Taxon object has the following attributes:
        name (str): canonical name of the taxon (e. g. `Bacillus subtilis').
        path (tuple): a sequence of taxa representing ancestors of this taxon from the root to the self. path[-1] is self, path[-2] is the parent.
        qps (boolean): either this taxon has QPS status.
    '''
    def __init__(self, identifier, base_url=None):
        '''Creates a Taxon object.

        Parameters:
            identifier (str): object identifier.
            base_url (str): base URL of the API. Defaults to None. When this parameter is None, then use this module BASE_URL value.
        '''
        LazyObject.__init__(self, identifier, base_url)

    def _object_url(self):
        '''/get/taxon/{self.identifier}'''
        return f'/get/taxon/{self.identifier}'

    def _fill(self, json):
        self.name = json['name']
        self.path = self.convert_path(json['path'][0])
        self.qps = (json['qps'] == 'yes')

    def __str__(self):
        return f'Taxon({self.identifier}, {self.name}, {LazyObject.path_string(self.path)}, QPS={self.qps})'

    @staticmethod
    def _identifier_param(taxon):
        if isinstance(taxon, Taxon):
            return taxon.identifier
        return str(taxon)

    @staticmethod
    def search(s=None, root=None, qps=None, base_url=None):
        '''Searches for taxa meeting criteria.

        Parameters:
            s (str): pattern to match taxon names or synonyms. POSIX-style regular expression.
            root (str|Taxon): if not None, only search for descendants of the specified taxon. Given as a taxon identifier or a Taxon object.
            qps (boolean): if True, only search for taxa that have QPS status.
            base_url (str): base URL of the API. Defaults to None. When this parameter is None, then use this module BASE_URL value.

        Returns:
            a generator of Taxon objects that match the query.
        '''
        params = {}
        if s is not None:
            params['s'] = s
        if root is not None:
            params['root'] = Taxon._identifier_param(root)
        if qps:
            params['qps'] = 'yes'
        for json in fetch_json('/search/taxon', base_url=base_url, params=params):
            yield Taxon(json['id'])


class OBTType(object):
    '''Types of OntoBiotope concepts. There are only 3 instances: HABITAT, PHENOTYPE, and USE.'''
    ALL = []
    MAP = {}

    def __init__(self, name):
        self.name = name
        if name in OBTType.MAP:
            raise ValueError('duplicate OBTType ' + name)
        OBTType.MAP[name] = self
        OBTType.ALL.append(self)

    @staticmethod
    def _name_param(obt_type):
        if isinstance(obt_type, OBTType):
            return obt_type.name
        return str(obt_type)


HABITAT = OBTType('habitat')
PHENOTYPE = OBTType('phenotype')
USE = OBTType('use')


class OBT(LazyObject):
    '''A OBT object represents a concept of the OntoBiotope ontology in the Omnicrobe database.

    Attributes:
        Inherits from LazyObject.
        When the data is loaded from the API, then the Taxon object has the following attributes:
        name (str): canonical name of the concept (e. g. `food').
        synonyms (tuple): other synonyms for this concept.
        paths (tuple): a collection of sequence of concepts representing ancestors of this concept from the root to the self.
        obt_type (OBTType): the concept type (habitat, phenotype or use).
    '''
    def __init__(self, identifier, base_url=None):
        LazyObject.__init__(self, identifier, base_url)

    def _object_url(self):
        '''/get/obt/{self.identifier}'''
        return f'/get/obt/{self.identifier}'

    def _fill(self, json):
        self.name = json['name']
        self.paths = tuple(self.convert_path(pathstr) for pathstr in json['path'])
        self.synonyms = json['synonyms'] if 'synonyms' in json else ()
        self.obt_type = OBTType.MAP[json['type']]

    def __str__(self):
        return f'OBT({self.identifier}, {self.obt_type.name}, {self.name}, {self.synonyms}, {tuple(LazyObject.path_string(p) for p in self.paths)})'

    @staticmethod
    def _identifier_param(obt):
        if isinstance(obt, OBT):
            return obt.identifier
        return str(obt)

    @staticmethod
    def search(s=None, root=None, obt_type=None, base_url=None):
        '''Searches for OntoBiotope concepts meeting criteria.

        Parameters:
            s (str): pattern to match concept names or synonyms. POSIX-style regular expression.
            root (str|OBT): if not None, only search for descendants of the specified concept. Given as a concept identifier or a OBT object.
            obt_type (str|OBTType): if not None, only search for concepts of the specified type. Given as a type name (habitat, phenotype or use) or a OBTType object.
            base_url (str): base URL of the API. Defaults to None. When this parameter is None, then use this module BASE_URL value.

        Returns:
            A generator of OBT objects that match the query.
        '''
        params = {}
        if s is not None:
            params['s'] = s
        if root is not None:
            params['root'] = OBT._identifier_param(root)
        if obt_type is not None:
            params['type'] = OBTType._name_param(obt_type)
        for json in fetch_json('/search/obt', base_url=base_url, params=params):
            yield OBT(json['id'])


class Relation:
    '''A Relation object represents a relation in Omnicrobe between a taxon and an OntoBiotope concept.

    Attributes:
        taxon (Taxon): the taxon side of the relation.
        obt (OBT): the OntoBiotope side of the relation.
        obt_type (OBTType): type of the concept on the OntoBiotope side of the relation. Is equal to self.obt.obt_type.
        source (str): name of the source where the relation was extracted.
        docs (list): list of document identifiers where the relation was extracted.
        taxon_forms (list): list of surface forms of the taxon.
        obt_forms (list): list of surface forms of the OntoBiotope concept.
    '''
    def __init__(self, taxon, obt, obt_type, source, docs, taxon_forms, obt_forms):
        self.taxon = taxon
        self.obt = obt
        self.obt_type = obt_type
        self.source = source
        self.docs = docs
        self.taxon_forms = taxon_forms
        self.obt_forms = obt_forms

    def __str__(self):
        return f'Relation({self.taxon.identifier} "{self.taxon.name}", {self.obt.identifier} "{self.obt.name}", {self.obt_type.name}, {self.source}, {self.docs}, {self.taxon_forms}, {self.obt_forms})'

    @staticmethod
    def search(sources=None, taxon=None, qps=None, obt=None, obt_type=None, base_url=None):
        '''Searches relations in the Omnicrobe database.

        Parameters:
            sources (str|seq): if not None, name of the source(s) where to search for relations.
            taxon (str|Taxon): if not None, only search for relations where the taxon is a descendant of the specified taxon. Given as a taxon identifier or a Taxon object.
            qps (boolean): if True only search for relations where the taxon has QPS status.
            obt (str|OBT): if not None, only search for relations where the OntoBiotope is a descendant of the specified concept. Given as a concept identifier or a OBT object.
            obt_type (str,OBTType): if not None, only search for relations where the OntoBiotope has the specified type. Given as a type name (habitat, phenotype or use) or a OBTType object.
            base_url (str): base URL of the API. Defaults to None. When this parameter is None, then use this module BASE_URL value.
        '''
        params = {}
        if sources is not None:
            if isinstance(sources, str):
                params['source'] = sources
            else:
                try:
                    params['source'] = ','.join(sources)
                except TypeError:
                    params['source'] = str(sources)
        if taxon is not None:
            params['taxid'] = Taxon._identifier_param(taxon)
        if qps:
            params['qps'] = 'yes'
        if obt is not None:
            params['obtid'] = OBT._identifier_param(obt)
        if obt_type is not None:
            try:
                params['type'] = ','.join(OBTType._name_param(t) for t in obt_type)
            except TypeError:
                params['type'] = OBTType._name_param(obt_type)
        for json in fetch_json('/search/relations', base_url=base_url, params=params):
            yield Relation(
                taxon=Taxon(json['taxid']),
                obt=OBT(json['obtid']),
                obt_type=OBTType.MAP[json['type']],
                source=json['source'],
                docs=json['docs'],
                taxon_forms=json['taxon_forms'],
                obt_forms=json['obt_forms']
            )


if __name__ == '__main__':
    print('Version')
    print(version())
    print('\n')

    print('taxon info by id')
    bs = Taxon('ncbi:1423')
    print(bs)
    print('\n')

    print('obt info by id')
    soil = OBT('OBT:000427')
    print(soil)
    print('\n')

    print('search taxon descendants')
    for t in Taxon.search(root=bs):
        print(t.identifier)
    print('\n')

    print('search taxon name')
    for t in Taxon.search('propionibacterium'):
        print(t.identifier)
    print('\n')

    print('search obt descendants')
    for o in OBT.search(root=soil):
        print(o.identifier)
    print('\n')

    print('search obt name')
    for o in OBT.search('food'):
        print(o.identifier)
    print('\n')

    print('search relations')
    for r in Relation.search(taxon=bs, obt=soil, qps=True, sources='GenBank'):
        print(r)
