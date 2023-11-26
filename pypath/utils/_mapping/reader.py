#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from __future__ import annotations

import pypath.share.session as _session


class MapReader(_session.Logger):
    """
    Reads ID translation data and creates ``MappingTable`` instances.
    When initializing ID conversion tables for the first time
    data is downloaded from UniProt and read into dictionaries.
    It takes a couple of seconds. Data is saved to pickle
    dumps, this way later the tables load much faster.
    """

    def __init__(
            self,
            param,
            ncbi_tax_id = None,
            entity_type = None,
            load_a_to_b = True,
            load_b_to_a = False,
            uniprots = None,
            lifetime = 300,
            resource_id_types = None,
        ):
        """
        Args
            param (MappingInput): A mapping table definition, any child of
                the `internals.input_formats.MappingInput` class.
            ncbi_tax_id (int): NCBI Taxonomy identifier of the organism.
            entity_type (str): An optional, custom string showing the type of
                the entities,  e.g. `protein`. This is not mandatory for the
                identification of mapping tables, hence the same name types
                can't be used for different entities. E.g. if both proteins
                and miRNAs have Entrez gene IDs then these should be
                different ID types (e.g. `entrez_protein` and `entrez_mirna`)
                or both protein and miRNA IDs can be loaded into one mapping
                table and simply called `entrez`.
            load_a_to_b (bool): Load the mapping table for translation from
                `id_type` to `target_id_type`.
            load_b_to_a (bool): Load the mapping table for translation from
                `target_id_type` to `id_type`.
            uniprots (set): UniProt IDs to query in case the source of the
                mapping table is the UniProt web service.
            lifetime (int): If this table has not been used for longer than
                this preiod it is to be removed at next cleanup. Time in
                seconds. Passed to ``MappingTable``.
            resource_id_types: Additional mappings between pypath and resource
                specific identifier type labels.
        """

        session_mod.Logger.__init__(self, name = 'mapping')

        self.ncbi_tax_id = (
            ncbi_tax_id or
            param.ncbi_tax_id or
            settings.get('default_organism')
        )

        self._log(
            'Reader created for ID translation table, parameters: '
            '`ncbi_tax_id=%u, id_a=%s, id_b=%s, '
            'load_a_to_b=%u, load_b_to_a=%u, '
            'input_type=%s (%s)`.' % (
                self.ncbi_tax_id,
                param.id_type_a,
                param.id_type_b,
                load_a_to_b,
                load_b_to_a,
                param.type,
                param.__class__.__name__,
            )
        )

        self.cachedir = cache_mod.get_cachedir()

        self.id_type_a = param.id_type_a
        self.id_type_b = param.id_type_b
        self.load_a_to_b = load_a_to_b
        self.load_b_to_a = load_b_to_a
        self.entity_type = entity_type
        self.source_type = param.type
        self.param = param
        self.lifetime = lifetime
        self.a_to_b = None
        self.b_to_a = None
        self.uniprots = uniprots
        self._resource_id_types = resource_id_types

        self.load()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def load(self):
        """
        The complete process of loading mapping tables. First sets up the
        paths of the cache files, then loads the tables from the cache files
        or the original sources if necessary. Upon successful loading from an
        original source writes the results to cache files.
        """

        self.use_cache = settings.get('mapping_use_cache')
        self.setup_cache()

        if self.use_cache:

            self.read_cache()

        if not self.tables_loaded():

            # read from the original source
            self.read()

            if self.tables_loaded():

                # write cache only at successful loading
                self.write_cache()


    @property
    def mapping_table_a_to_b(self):
        """
        Returns a ``MappingTable`` instance created from the already
        loaded data.
        """

        return self._get_mapping_table('a', 'b')


    @property
    def mapping_table_b_to_a(self):
        """
        Returns a ``MappingTable`` instance created from the already
        loaded data.
        """

        return self._get_mapping_table('b', 'a')


    def id_type_side(self, id_type):
        """
        Tells if an ID type is on the "a" or "b" (source or target) side
        in the current mapping table definition.

        Args
            id_type (str): An ID type label.

        Returns
            Returns the string "a" if `id_type` is on the source side in
            the mapping table definition, "b" if it is on the target side,
            None if the `id_type` is not in the definition.
        """

        return (
            'a'
                if id_type == self.id_type_a else
            'b'
                if id_type == self.id_type_b else
            None
        )


    def _get_mapping_table(self, *args):

        data = getattr(self, '%s_to_%s' % args)
        id_type = getattr(self, 'id_type_%s' % args[0])
        target_id_type = getattr(self, 'id_type_%s' % args[1])

        if isinstance(data, dict):

            return MappingTable(
                data = data,
                id_type = id_type,
                target_id_type = target_id_type,
                ncbi_tax_id = self.ncbi_tax_id,
                lifetime = self.lifetime,
            )


    def tables_loaded(self):
        """
        Tells if the requested tables have been created.
        """

        return (
            (bool(self.a_to_b) or not self.load_a_to_b) and
            (bool(self.b_to_a) or not self.load_b_to_a)
        )


    def write_cache(self):
        """
        Exports the ID translation data into pickle files.
        """

        self._write_cache('a', 'b')
        self._write_cache('b', 'a')


    def _write_cache(self, *args):

        data = getattr(self, '%s_to_%s' % args)

        if self._to_be_loaded(*args) and data:

            cachefile = self._attr('cachefile', *args)

            self._remove_cache_file(*args)

            pickle.dump(data, open(cachefile, 'wb'))


    def read_cache(self):
        """
        Reads the ID translation data from a previously saved pickle file.
        """

        self._read_cache('a', 'b')
        self._read_cache('b', 'a')


    def _read_cache(self, *args):

        if self._to_be_loaded(*args):

            cachefile = self._attr('cachefile', *args)

            if os.path.exists(cachefile):

                with open(cachefile, 'rb') as fp:

                    from_cache = pickle.load(fp)

                setattr(
                    self,
                    '%s_to_%s' % args,
                    from_cache,
                )
                self._log(
                    'Loading `%s` to `%s` mapping table '
                    'from pickle file `%s`.' % (
                        self.param.id_type_a,
                        self.param.id_type_b,
                        cachefile,
                    )
                )


    def _to_be_loaded(self, *args):

        return self._attr('load', *args)


    def _attr(self, attr, *args):

        return getattr(self, self._attr_name(attr, *args))


    @staticmethod
    def _attr_name(attr, *args):

        return '%s_%s_to_%s' % ((attr,) + args)


    def read(self):
        """
        Reads the ID translation data from the original source.
        """

        method = 'read_mapping_%s' % self.source_type

        if hasattr(self, method):

            getattr(self, method)()


    def setup_cache(self):
        """
        Constructs the cache file path as md5 hash of the parameters.
        """

        self._setup_cache('a', 'b')
        self._setup_cache('b', 'a')


    def _setup_cache(self, *args):

        mapping_id_attr = self._attr_name('mapping_id', *args)
        cachefile_attr  = self._attr_name('cachefile', *args)

        setattr(
            self,
            mapping_id_attr,
            self._get_mapping_id(*args),
        )

        setattr(
            self,
            cachefile_attr,
            os.path.join(self.cachedir, getattr(self, mapping_id_attr)),
        )


    def _get_mapping_id(self, *args):
        """
        Returns an md5 checksum unambigously identifying the mapping table
        by the identifiers, the direction of translation, the organism
        and other parameters like, for example, the source URL.
        """

        return common.md5(
            json.dumps(
                (
                    getattr(self, 'id_type_%s' % args[0]),
                    getattr(self, 'id_type_%s' % args[1]),
                    self.ncbi_tax_id,
                    sorted(self.param.__dict__.items())
                )
            )
        )


    def _cache_files_exist(self):
        """
        Checks if both cache files are either not necessary or exist.
        """

        return (
            self.cache_file_exists('a', 'b') and
            self.cache_file_exists('b', 'a')
        )


    def _cache_file_exists(self, *args):
        """
        Checks if a cache file is either not necessary or exists.
        """

        return (
            not self._attr('load', *args) or
            os.path.isfile(self._attr('cachefile', *args))
        )


    def _remove_cache_file(self, *args):

        cachefile = self._attr('cachefile', *args)

        if os.path.exists(cachefile):

            self._log('Removing mapping table cache file `%s`.' % cachefile)
            os.remove(cachefile)


    def read_mapping_file(self):
        """
        Reads a mapping table from a local file or a function.
        """

        if not os.path.exists(self.param.input):

            method = inputs.get_method(self.param.input)

            if not method:

                return {}

            else:

                input_args = (
                    self.param.input_args
                        if hasattr(self.param, 'input_args') else
                    {}
                )
                infile = method(**input_args)

        else:

            infile = open(self.param.input, encoding = 'utf-8', mode = 'r')
            total = os.path.getsize(self.param.input)

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)

        for i, line in enumerate(infile):

            if self.param.header and i < self.param.header:

                continue

            if hasattr(line, 'decode'):

                line = line.decode('utf-8')

            if hasattr(line, 'rstrip'):

                line = line.rstrip().split(self.param.separator)

            if len(line) < max(self.param.col_a, self.param.col_b):

                continue

            id_a = line[self.param.col_a]
            id_b = line[self.param.col_b]

            if self.load_a_to_b:

                a_to_b[id_a].add(id_b)

            if self.load_b_to_a:

                b_to_a[id_b].add(id_a)

        if hasattr(infile, 'close'):

            infile.close()

        self.a_to_b = a_to_b if self.load_a_to_b else None
        self.b_to_a = b_to_a if self.load_b_to_a else None


    @staticmethod
    def _uniprotkb_id_type(id_type: str) -> bool:

        return input_formats.UniprotListMapping._uniprotkb_id_type(
            id_type,
        )


    def read_mapping_uniprot_list(self):
        """
        Builds a mapping table by downloading data from UniProt's
        upload lists service.
        """

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)
        swap = False

        if not self.uniprots:

            self.set_uniprot_space()

        # We need a list to query this service, and we have method only for
        # getting a proteome wide list of UniProt IDs. If the translated
        # ID type is not UniProt, then first we need to translate the
        # proteome wide reference list from UniProt to the target ID type.
        if not self._uniprotkb_id_type(self.param.id_type_a):

            if self._uniprotkb_id_type(self.param.id_type_b):

                swap = True
                self.param.swap_sides()
                self.load_a_to_b, self.load_b_to_a = (
                    self.load_b_to_a,
                    self.load_a_to_b,
                )
                upload_ac_list = self.uniprots

            else:

                u_target = self._read_mapping_uniprot_list(
                    uniprot_id_type_a = 'UniProtKB_AC-ID',
                    uniprot_id_type_b = self.param.uniprot_id_type_a,
                )

                upload_ac_list = [l.split('\t')[1].strip() for l in u_target]

        else:

            upload_ac_list = self.uniprots

        uniprot_data = self._read_mapping_uniprot_list(
            upload_ac_list = upload_ac_list,
        )
        ens = (
            self.param.id_type_a.startswith('ens') or
            self.param.id_type_b.startswith('ens') or
            'ensembl' in self.param.id_type_a.lower() or
            'ensembl' in self.param.id_type_b.lower()
        )
        reens = re.compile(r'(ENS[A-Z]+\d+)\.\d+')

        for l in uniprot_data:

            if not l:

                continue

            if ens:

                l = reens.sub(r'\1', l)

            l = l.strip().split('\t')

            if self.load_a_to_b:

                a_to_b[l[0]].add(l[1])

            if self.load_b_to_a:

                b_to_a[l[1]].add(l[0])

        if swap:

            a_to_b, b_to_a = b_to_a, a_to_b
            self.load_a_to_b, self.load_b_to_a = (
                self.load_b_to_a,
                self.load_a_to_b,
            )
            self.param.swap_sides()

        self.a_to_b = a_to_b if self.load_a_to_b else None
        self.b_to_a = b_to_a if self.load_b_to_a else None


    def set_uniprot_space(self, swissprot = None):
        """
        Sets up a search space of UniProt IDs.

        Args
            swissprot (bool): Use only SwissProt IDs, not TrEMBL. True
                loads only SwissProt IDs, False only TrEMBL IDs, None
                loads both.
        """

        swissprot = self.param.swissprot if swissprot is None else swissprot

        self.uniprots = uniprot_db.all_uniprots(
            self.ncbi_tax_id,
            swissprot = swissprot,
        )


    def _read_mapping_uniprot_list(
            self,
            uniprot_id_type_a = None,
            uniprot_id_type_b = None,
            upload_ac_list = None,
            chunk_size = None,
        ):
        """
        Reads a mapping table from UniProt "upload lists" service.

        Args
            uniprot_id_type_a (str): Source ID type label as used in UniProt.
            uniprot_id_type_b (str): Target ID type label as used in UniProt.
            upload_ac_list (list): The identifiers to use in the query to
                the ID Mapping service. By default the list of all UniProt
                IDs for the organism is used.
            chunk_size (int): Number of IDs in one query. Too large queries
                might fail, by default we include 100,000 IDs in one query.
        """

        chunk_size = (
            chunk_size or
            settings.get('uniprot_uploadlists_chunk_size')
        )
        uniprot_id_type_a = uniprot_id_type_a or self.param.uniprot_id_type_a
        uniprot_id_type_b = uniprot_id_type_b or self.param.uniprot_id_type_b

        if not upload_ac_list:

            self._log(
                'No identifiers provided, '
                'using all UniProt IDs of the organism.'
            )
            upload_ac_list = self.uniprots

        upload_ac_list = sorted(upload_ac_list)

        self._log(
            'Querying the UniProt ID Mapping service for ID translation '
            'data. Querying a list of %u IDs.' % len(upload_ac_list)
        )

        run_url = urls.urls['uniprot_idmapping']['run']
        poll_result = {}
        result = []

        # loading data in chunks of 10,000 by default
        for i in range(math.ceil(len(upload_ac_list) / chunk_size)):

            this_chunk = upload_ac_list[i * chunk_size:(i + 1) * chunk_size]

            self._log(
                'Request to UniProt ID Mapping, chunk #%u with %u IDs.' % (
                    i,
                    len(this_chunk),
                )
            )

            post = {
                'from': uniprot_id_type_a,
                'to': uniprot_id_type_b,
                'ids': ' '.join(sorted(this_chunk)),
            }
            accept_json = {'req_headers': ['Accept: application/json']}

            run_args = {'url': run_url, 'post': post}
            nocache = {'cache': False, 'large': False}
            large = {'silent': False, 'large': True}

            cache_path = curl.Curl.cache_path(**run_args)

            if not os.path.exists(cache_path):

                run_c = curl.Curl( **run_args, **nocache, **accept_json)

                if run_c.status != 200:

                    raise RuntimeError(
                        'Failed to submit job to UniProt ID Mapping. '
                        'See details in the log.'
                    )

                jobid = json.loads(run_c.result)['jobId']

                self._log(
                    f'Submitted job to UniProt ID Mapping, job ID: `{jobid}`.'
                )

                timeout = settings.get('uniprot_idmapping_timeout')
                interval = settings.get('uniprot_idmapping_poll_interval')
                max_polls = math.ceil(timeout / interval)
                poll_url = urls.urls['uniprot_idmapping']['poll'] % jobid
                poll_args = {'url': poll_url} | nocache | accept_json

                for i in range(max_polls):

                    self._log(
                        f'Polling job UniProt ID Mapping job `{jobid}`, '
                        f'poll {i + 1} of {max_polls}.'
                    )

                    poll_c = curl.Curl(**poll_args)

                    if poll_c.status != 200:

                        self._log(f'Poll failed with HTTP {poll_c.status}.')
                        continue

                    poll_result = json.loads(poll_c.result)

                    if 'status' in poll_result or 'failedIds' in poll_result:

                        self._log(
                            f'UniProt ID Mapping job `{jobid}` '
                            'successfully completed.'
                        )
                        break

                    elif 'messages' in poll_result:

                        msg = (
                            'UniProt ID Mapping job failed: ' +
                            ' '.join(common.to_list(poll_result['messages']))
                        )

                        self._log(msg)

                        raise RuntimeError(msg)

                    time.sleep(interval)

                self._log(
                    'Getting UniProt ID Mapping results URL '
                    'for job `{jobid}`.'
                )
                det_url = urls.urls['uniprot_idmapping']['details'] % jobid
                det_c = curl.Curl(url = det_url, **nocache, **accept_json)
                result_url = (
                    json.loads(det_c.result)['redirectURL'].
                    replace('/idmapping/results/', '/idmapping/stream/').
                    replace('/results/', '/results/stream/').
                    __add__('?format=tsv')
                )

                self._log(
                    'Retrieving UniProt ID Mapping results '
                    f'from `{result_url}`.'
                )

                with curl.cache_delete_on():

                    res_c = curl.Curl(
                        url = result_url,
                        cache = cache_path,
                        **large
                    )

            else:

                res_c = curl.Curl(**run_args, **large)

            result.extend(list(res_c.fileobj)[1:])

        return result


    def read_mapping_uniprot(self):
        """
        Downloads ID mappings directly from UniProt.
        See the names of possible identifiers here:
        http://www.uniprot.org/help/programmatic_access
        """

        query = uniprot_input.UniprotQuery(
            reviewed = True if self.param.swissprot else None,
            organism = self.ncbi_tax_id,
            fields = self.param._resource_id_type_a,
        )
        self._log(f'UniProt REST API call: `{query.url_plain}`.')
        trembl = 'trembl' in self.param
        protein_name = self.param.field == 'protein names'
        query.name_process = not protein_name and not trembl
        data = query.perform()

        if not query.name_process:

            def maybe_split(v):

                if trembl and not any(ch.islower() for ch in v):
                    v = common.del_empty(query._FIELDSEP.split(v))
                elif protein_name:
                    v = self._process_protein_name(v)

                return v


            data = {k: maybe_split(v) for k, v in data.items()}

        data = {k: common.to_set(v) for k, v in data.items()}

        self.a_to_b = (
            common.swap_dict(data, force_sets = True)
                if self.load_a_to_b else
            None
        )
        self.b_to_a = data if self.load_b_to_a else None


    def read_mapping_pro(self):

        pro_data = pro_input.pro_mapping(target_id_type = self.param.id_type)

        pro_to_other = collections.defaultdict(set)

        for pro, other in pro_data:

            pro_to_other[pro].add(other)

        self.a_to_b = (
            None
                if not self.load_a_to_b else
            common.swap_dict(pro_to_other, force_sets = True)
                if self.param.to_pro else
            dict(pro_to_other)
        )
        self.b_to_a = (
            None
                if not self.load_b_to_a else
            dict(pro_to_other)
                if self.param.to_pro else
            common.swap_dict(pro_to_other, force_sets = True)
        )


    def read_mapping_biomart(self):
        """
        Loads a mapping table using BioMart data.
        """

        ens_organism = taxonomy.ensure_ensembl_name(self.param.ncbi_tax_id)

        if not ens_organism:

            self._log(
                'Organism not available in Ensembl: `%u`.' % (
                    self.param.ncbi_tax_id
                )
            )
            return

        dataset = '%s_gene_ensembl' % ens_organism
        biomart_data = biomart_input.biomart_query(
            attrs = self.param.attrs,
            dataset = dataset,
        )

        a_to_b = collections.defaultdict(set)
        b_to_a = collections.defaultdict(set)

        for rec in biomart_data:

            id_a = getattr(rec, self.param.biomart_id_type_a)
            id_b = getattr(rec, self.param.biomart_id_type_b)

            if id_a and id_b:

                if self.load_a_to_b:

                    a_to_b[id_a].add(id_b)

                if self.load_b_to_a:

                    b_to_a[id_b].add(id_a)

        self.a_to_b = dict(a_to_b) if self.load_a_to_b else None
        self.b_to_a = dict(b_to_a) if self.load_b_to_a else None


    def read_mapping_array(self):
        """
        Loads mapping table between microarray probe IDs and genes.
        """

        probe_mapping = biomart_input.biomart_microarrays(
            organism = self.param.ncbi_tax_id,
            vendor = self.param.array_id,
            gene = self.param.ensembl_id == 'ensg',
            transcript = self.param.ensembl_id == 'enst',
            peptide = self.param.ensembl_id == 'ensp',
        )

        a_to_b__probe_to_gene = self.param.id_type_a == self.param.array_id

        if (
            (
                a_to_b__probe_to_gene and
                self.load_a_to_b
            ) or (
                not a_to_b__probe_to_gene and
                self.load_b_to_a
            )
        ):

            probe_to_gene = collections.defaultdict(set)

            for ensembl_id, probes in iteritems(probe_mapping):

                for probe in probes:

                    probe_to_gene[probe.probe].add(ensembl_id)

            setattr(
                self,
                'a_to_b' if a_to_b__probe_to_gene else 'b_to_a',
                dict(probe_to_gene),
            )

        if (
            (
                a_to_b__probe_to_gene and
                self.load_b_to_a
            ) or (
                not a_to_b__probe_to_gene and
                self.load_a_to_b
            )
        ):

            gene_to_probe = dict(
                (
                    ensembl_id,
                    {p.probe for p in probe_ids}
                )
                for ensembl_id, probe_ids in iteritems(probe_mapping)
            )

            setattr(
                self,
                'b_to_a' if a_to_b__probe_to_gene else 'a_to_b',
                gene_to_probe,
            )


    def _read_mapping_smallmolecule(self):
        """
        Loads a small molecule ID translation table.
        """

        if self.param.input_method:

            method = inputs.get_method(self.param.input_method)

        else:

            mod = globals()[f'{self.source_type}_input']
            method = getattr(mod, f'{self.source_type}_mapping')

        data = method(
            id_type_a = self.resource_id_type_a,
            id_type_b = self.resource_id_type_b,
        )

        if self.load_a_to_b:

            self.a_to_b = data

        if self.load_b_to_a:

            self.b_to_a = common.swap_dict(data, force_sets = True)

        self.ncbi_tax_id = _const.NOT_ORGANISM_SPECIFIC


    def read_mapping_ramp(self):
        """
        Loads an ID translation table from RaMP.
        """

        self._read_mapping_smallmolecule()


    def read_mapping_unichem(self):
        """
        Loads an ID translation table from UniChem.
        """

        self._read_mapping_smallmolecule()


    def read_mapping_hmdb(self):
        """
        Loads an ID translation table from th Human Metabolome Database.
        """

        self._read_mapping_smallmolecule()


    @staticmethod
    def _process_protein_name(name):

        rebr = re.compile(r'\(([^\)]{3,})\)')
        resq = re.compile(r'\[([^\]]{3,})\]')

        names = [name.split('(')[0]]
        names += rebr.findall(name)
        others = common.flat_list([x.split(';') for x in resq.findall(name)])
        others = [x.split(':')[1] if ':' in x else x for x in others]
        others = [x.split('(')[1] if '(' in x else x for x in others]
        names += others

        return {x.strip() for x in names}


    def resource_id_type(self, side = Literal['a', 'b']) -> str | None:
        """
        Resource specific identifier type.
        """

        return (
            getattr(self.param, f'resource_id_type_{side}') or
            self._resource_id_types.get(getattr(self.param, f'id_type_{side}'))
        )


    @property
    def resource_id_type_a(self) -> str | None:

        return self.resource_id_type('a')


    @property
    def resource_id_type_b(self) -> str | None:

        return self.resource_id_type('b')
