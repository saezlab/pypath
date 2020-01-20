
# How to load TFregulons in `pypath`

TFregulons is a comprehensive collection of transcription factor regulons combining evidences of four main categories: 1) literature curation, 2) ChIP-Seq, 3) PWM based in silico TFBS matching and 4) reverse engineering from expression data. TFregulons provide a confidence score for the interactions based on TF activity prediction benchmarks. The levels of the score go from `A` (supported by most reliable evidences) to `E` (supported by the weakest level of evidences). `pypath` is able to load any combination of the levels. Note, level `D` and `E` contain already a huge number of interactions which require high amount of memory to load. In addition it is possible to limit interactions only to the literature curated ones using the `only_curated = True` argument. Literature references are not available from TFregulons. For more details see Garcia-Alonso 2018 and https://github.com/saezlab/DoRothEA.

In `pypath` the `data_formats.transcription` input definition set now contains only TFregulons as this is already the most comprehensive integrated resource of transcriptional regulation.


```python
import pypath
```


```python
pypath.data_formats.transcription
```




    {'tfregulons': <pypath.input_formats.ReadSettings at 0x688bf7fd8c88>}



However many of the literature curated databases are covered also in the `data_formats.transcription_onebyone` input definition set.


```python
pypath.data_formats.transcription_onebyone
```




    {'abs': <pypath.input_formats.ReadSettings at 0x688bf7fd8b00>,
     'encode_dist': <pypath.input_formats.ReadSettings at 0x688bf7fd8b38>,
     'encode_prox': <pypath.input_formats.ReadSettings at 0x688bf7fd8b70>,
     'htri': <pypath.input_formats.ReadSettings at 0x688bf7fd8be0>,
     'oreganno': <pypath.input_formats.ReadSettings at 0x688bf7fd8c18>,
     'pazar': <pypath.input_formats.ReadSettings at 0x688bf7fd8ba8>,
     'signor': <pypath.input_formats.ReadSettings at 0x688bf7fd8c50>}



By default TF-target relationships with `A` and `B` confidence levels loaded.


```python
pa = pypath.PyPath()
```

    
    	=== d i s c l a i m e r ===
    
    	All data coming with this module
    	either as redistributed copy or downloaded using the
    	programmatic interfaces included in the present module
    	are available under public domain, are free to use at
    	least for academic research or education purposes.
    	Please be aware of the licences of all the datasets
    	you use in your analysis, and please give appropriate
    	credits for the original sources when you publish your
    	results. To find out more about data sources please
    	look at `pypath.descriptions` and
    	`pypath.data_formats.urls`.
    
    	> New session started,
    	session ID: 'gol2f'
    	logfile: './log/gol2f.log'
    	pypath version: 0.7.71



```python
pa.load_tfregulons()
```

    	:: Loading data from cache previously downloaded from www.uniprot.org
    	:: Ready. Resulted `plain text` of type unicode string.                                                                                              
    	:: Local file at `/home/denes/Dokumentumok/pw/dev/update_cache/cache/ec920965677ac83b8805d72853c79d45-`.
     > TFRegulons
    	:: Loading data from cache previously downloaded from saezlab.org
    	:: Ready. Resulted `plain text` of type file object.                                                                                                 
    	:: Local file at `/home/denes/Dokumentumok/pw/dev/update_cache/cache/0596f5303ebb6d4f4f9a81fc8110453f-tfregulons_database_v01_20180216__AB.tsv`.
    	:: Loading 'genesymbol' to 'uniprot' mapping table
    	:: Loading 'uniprot-sec' to 'uniprot-pri' mapping table
    	:: Loading data from cache previously downloaded from www.uniprot.org
    	:: Ready. Resulted `plain text` of type unicode string.                                                                                              
    	:: Local file at `/home/denes/Dokumentumok/pw/dev/update_cache/cache/842e6f2bc63115660aec8aff917330ce-`.
    	:: Loading data from cache previously downloaded from ftp.uniprot.org
    	:: Ready. Resulted `plain text` of type file object.                                                                                                 
    	:: Local file at `/home/denes/Dokumentumok/pw/dev/update_cache/cache/49314fe217bf0f2a5544a2c4314b4adf-sec_ac.txt`.


            Reading from file -- finished: 0.00it [00:00, ?it/s]


    	:: Loading 'genesymbol' to 'trembl' mapping table
    	:: Loading 'genesymbol-syn' to 'uniprot' mapping table


            Processing nodes -- finished: 100%|██████████| 14.3K/14.3K [00:00<00:00, 514Kit/s]
            Processing edges -- finished: 100%|██████████| 14.3K/14.3K [00:00<00:00, 168Kit/s]
            Processing attributes -- finished: 100%|██████████| 14.3K/14.3K [00:05<00:00, 2.79Kit/s]


    
     :: Comparing with reference lists... done.
    
     > 14080 interactions between 4877 nodes
     from 17 resources have been loaded,
     for details see the log: ./log/gol2f.log


You can set other confidence levels and also get only the literature curated interactions at any confidence level.


```python
pa = pypath.PyPath()
pa.load_tfregulons(levels = {'A', 'B', 'C'}, curated_only = True)
```

The type of these interactions is `TF`:


```python
pa.graph.es[0]['type']
```




    ['TF']




```python
pa.graph.es[0]['sources']
```




    {'HTRIdb', 'trrust_signed'}




```python
pa.graph.es[0]['tfregulons_curated']
```




    True



You can also access the TFregulons data as a `list` using the `pypath.dataio.get_tfregulons` method.


```python
tfregulons = list(pypath.dataio.get_tfregulons(levels = {'A'}))
```

    	:: Loading data from cache previously downloaded from saezlab.org
    	:: Ready. Resulted `plain text` of type file object.                                                                                                 
    	:: Local file at `/home/denes/Dokumentumok/pw/dev/update_cache/cache/85ca87819310f9f39fefda4327748bc9-tfregulons_database_v01_20180216__A.tsv`.



```python
tfregulons[0]
```




    ['AHR',
     'CYP1A1',
     '0',
     'A',
     True,
     False,
     True,
     False,
     'HTRIdb,trrust_signed',
     '',
     'hocomoco_v11',
     '',
     'HTRIdb,trrust_signed,hocomoco_v11']


