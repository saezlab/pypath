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

import pypath.share.common as common
import pypath.internals.annot_formats as af


#TODO should go to jsons
"""
Gene Ontology annotations to select categories relevant in intercellular
signaling.
"""
go_combined_classes = {
    'junction':
        """
        cell junction OR
        cell junction assembly OR
        cell junction organization OR
        intercellular bridge organization OR
        gap junction-mediated intercellular transport OR
        gap junction
        """,
    'gap junction':
        """
        gap junction
        """,
    'tight junction':
        """
        tight junction
        """,
    'cell-substrate junction':
        """
        cell-substrate junction
        """,
    'cell-cell junction':
        """
        cell-cell junction
        """,
    'extracellular':
        """
        extracellular region
        """,
    'intracellular':
        """
        intracellular organelle OR
        intracellular organelle lumen OR
        intracellular anatomical structure
        """,
    'cell_surface':
        """
        cell surface OR
        external side of plasma membrane
        """,
    'transmembrane':
        """
        integral component of membrane OR
        transmembrane signaling receptor activity
        """,
    'ecm':
        """
        extracellular matrix OR
        complex of collagen trimers OR
        collagen network OR
        banded collagen fibril OR
        collagen beaded filament OR
        elastic fiber OR
        fibronectin fibril
        """,
    'extracell enzyme':
        """
        catalytic activity AND
        extracellular region
        """,
    'enzyme':
        """
        catalytic activity
        """,
    'extracell peptidase':
        """
        peptidase activity AND
        extracellular region
        """,
    'peptidase':
        """
        peptidase activity
        """,
    'growth factor binding':
        """
        growth factor binding AND
        extracellular space
        """,
    'receptor regulation':
        """
        (signaling receptor regulator activity OR
        regulation of receptor recycling OR
        receptor clustering OR
        receptor diffusion trapping OR
        membrane raft localization) AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane OR
        intrinsic component of plasma membrane)
        """,
    'receptor inhibition':
        """
        signaling receptor inhibitor activity AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane OR
        intrinsic component of plasma membrane)
        """,
    'receptor activation':
        """
        signaling receptor activator activity AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane OR
        intrinsic component of plasma membrane)
        """,
    'ligands':
        """
        receptor ligand activity AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane)
        """,
    'secreted ligands':
        """
        receptor ligand activity AND
        extracellular region
        """,
    'surface ligands':
        """
        receptor ligand activity AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'receptors':
        """
        (signaling receptor activity OR
        cell surface receptor signaling pathway OR
        transmembrane signaling receptor activity OR
        receptor complex OR
        (cellular response to stimulus AND signal transduction)) AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'membrane ligands':
        """
        receptor ligand activity AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'hormone receptors':
        """
        hormone binding AND
        (cell surface OR
        external side of plasma membrane)
        """,
    'ecm structure':
        """
        (extracellular region OR
        extracellular matrix) AND
        (structural molecule activity OR
        extracellular matrix structural constituent OR
        structural constituent of bone)
        """,
    'ecm production':
        """
        extracellular matrix assembly OR
        extracellular matrix organization OR
        extracellular matrix constituent secretion OR
        collagen metabolic process
        """,
    'endocytosis':
        """
        cargo adaptor activity
        """,
    'adhesion to matrix':
        """
        cell adhesion mediator activity OR
        extracellular matrix binding OR
        hydroxyapatite binding OR
        cell-substrate adhesion
        """,
    'response to adhesion':
        """
        cellular response to cell-matrix adhesion OR
        contact inhibition OR
        establishment or maintenance of cell polarity OR
        extracellular matrix-cell signaling
        """,
    'adhesion to other cells':
        """
        cell-cell adhesion
        """,
    'adhesion':
        """
        cell adhesion mediator activity OR
        extracellular matrix binding OR
        hydroxyapatite binding OR
        cell adhesion OR
        cell-cell adhesion OR
        cell-substrate adhesion OR
        cell adhesion molecule production OR
        cell-cell adhesion in response to extracellular stimulus OR
        cellular response to cell-matrix adhesion OR
        contact inhibition OR
        establishment or maintenance of cell polarity
        """,
    'cell-cell signaling':
        """
        cell-cell signaling
        """,
    'transport':
        """
        intrinsic component of plasma membrane AND
        transmembrane transporter activity
        """,
    'regulation of transport':
        """
        (cell surface OR
        external side of plasma membrane OR
        extracellular region) AND
        (regulation of transmembrane transporter activity OR
        channel regulator activity)
        """,
    'ion channels':
        """
        (ion channel activity OR
        ion transmembrane transporter activity) AND
        intrinsic component of plasma membrane
        """,
    'regulation of ion channels':
        """
        (cell surface OR
        external side of plasma membrane OR
        extracellular region) AND
        (ion channel regulator activity OR
        regulation of ion transmembrane transporter activity OR
        channel regulator activity)
        """,
    'autocrine signaling':
        """
        autocrine signaling
        """,
    'paracrine signaling':
        """
        paracrine signaling
        """,
    'signal release':
        """
        signal release OR
        hormone secretion OR
        hormone metabolic process OR
        cytokine production
        """,
    'secretion':
        """
        exocytic process OR
        secretion by cell
        """,
    'communication by exosomes':
        """
        extracellular vesicle biogenesis OR
        extracellular exosome assembly OR
        vesicle-mediated intercellular transport OR
        cell-cell signaling via exosome
        """,
    'membrane docking':
        """
        membrane docking AND
        (extracellular region OR
        cell surface OR
        external side of plasma membrane)
        """
}


go_single_terms = {

    # cellular component
    'C': {
        # junction
        'cell junction',

        # extracellular
        'extracellular region',

        # extracellular_matrix
        'extracellular matrix',
        'complex of collagen trimers',
        'collagen network',
        'banded collagen fibril',
        'collagen beaded filament',
        'elastic fiber',
        'fibronectin fibril',

        # exosome
        # could not find sub-term for their membrane or lumen
        'extracellular vesicle',

        # cell_surface
        # only plasma membrane components facing outside
        'cell surface',
        'external side of plasma membrane',

        # these contains also intracellular components
        'plasma membrane',
        'extrinsic component of plasma membrane',
        'intrinsic component of plasma membrane',
        'cytoplasmic side of plasma membrane',

        'immunological synapse',

        'clathrin-coated pit',
        'plasma membrane raft',

        'presynaptic membrane',
        'presynaptic endocytic zone',
        'presynaptic endocytic zone membrane',
        'intrinsic component of presynaptic membrane',
        'extrinsic component of presynaptic membrane',
        'presynaptic active zone membrane',

        'postsynaptic density membrane',
        'intrinsic component of postsynaptic density membrane',
        'extrinsic component of postsynaptic density membrane',

        'synaptic vesicle',
        'intrinsic component of synaptic vesicle membrane',
        'synaptic vesicle lumen',
        'synaptic vesicle membrane',
        'extracellular matrix of synaptic cleft',

        'axolemma',
        'neuron projection membrane',
        'neuronal cell body membrane',
        'photoreceptor inner segment membrane',
        'photoreceptor outer segment membrane',
        'stereocilium membrane',
        'stereocilia coupling link',

    },

    # molecular function
    'F': {

        # upper level terms which help to categorize
        # molecular functions in the intercellular signaling
        'molecular carrier activity',
        'cargo receptor activity',
        'binding',
        'regulation of binding',
        'positive regulation of binding',
        'negative regulation of binding',
        'protein folding chaperone',
        'antioxidant activity',

        # upper level generic terms for regulation
        'regulation of molecular function',
        'negative regulation of molecular function',
        'positive regulation of molecular function',

        # enzymes
        'catalytic activity',
        'regulation of catalytic activity',
        'positive regulation of catalytic activity',
        'negative regulation of catalytic activity',
        'enzyme regulator activity',
        'enzyme activator activity',
        'enzyme inhibitor activity',
        'catalytic activity, acting on a protein',
        # peptidases
        'peptidase activity',
        'peptidase regulator activity',
        'peptidase activator activity',
        'peptidase inhibitor activity',
        'regulation of peptidase activity',
        'peptidase activator activity',
        # receptors
        'signaling receptor activator activity',
        'signaling receptor inhibitor activity',
        'receptor ligand activity',
        'neurotransmitter receptor regulator activity',
        'signaling receptor activity',
        'negative regulation of signaling receptor activity',
        'positive regulation of signaling receptor activity',
        'signaling receptor activator activity',
        'signaling receptor inhibitor activity',
        'regulation of signaling receptor activity',
        'receptor complex',
        (
            'transforming growth factor beta receptor,'
            'cytoplasmic mediator activity'
        ),
        # other relevant binding activities
        'antigen binding',
        'hormone binding',
        'neurotransmitter binding',

        # endocytosis
        'cargo adaptor activity',

        # ECM, structural proteins
        'structural molecule activity',
        'extracellular matrix structural constituent',
        'structural constituent of bone',

        # adhesion to base membrane and ECM
        'cell adhesion mediator activity',
        'extracellular matrix binding',
        'hydroxyapatite binding',

        # transporters
        'transporter activity',
        'regulation of transporter activity',
        'positive regulation of transporter activity',
        'negative regulation of transporter activity',
        'transmembrane transporter activity',
        'drug transmembrane transporter activity',
        'regulation of transmembrane transporter activity',
        # channels
        'channel regulator activity',
        'channel activator activity',
        'channel inhibitor activity',
        # ion channels
        'ion channel inhibitor activity',
        'ion channel regulator activity',
        'ion transmembrane transporter activity',
        'regulation of ion transmembrane transporter activity',
        'positive regulation of ion transmembrane transporter activity',
        'negative regulation of ion transmembrane transporter activity',
    },

    # biological process
    'P': {

        # cell adhesion
        'cell adhesion',
        'cell-cell adhesion',
        'cell-substrate adhesion',
        'cell adhesion molecule production',
        'cell-cell adhesion in response to extracellular stimulus',
        'cellular response to cell-matrix adhesion',
        'contact inhibition',
        'establishment or maintenance of cell polarity'

        # cellular responses
        'myofibroblast cell apoptotic process',
        'fibroblast apoptotic process',
        'epithelial cell apoptotic process',

        # cell activation
        'cell activation',
        'endothelial cell activation',
        'fibroblast activation',
        'leukocyte activation',

        # cell communication
        'cell communication',
        'cell communication by chemical coupling',
        'cell communication by electrical coupling',
        'cell-cell signaling',
        'autocrine signaling',
        'paracrine signaling',
        'cell-cell signaling via exosome',
        'epithelial-mesenchymal cell signaling',
        'cellular response to extracellular stimulus',
        'synaptic signaling',
        'cell-cell recognition',

        # signal release
        'signal release',
        'hormone secretion',
        'exocytic process',
        'secretion by cell',
        'hormone metabolic process',
        'cytokine production',
        'cytokine secretion',

        # receptors
        'regulation of receptor recycling',
        'receptor clustering',
        'receptor diffusion trapping',
        'membrane raft localization',

        # exosomes
        'extracellular vesicle biogenesis',
        'extracellular exosome assembly',
        'vesicle-mediated intercellular transport',

        # junction
        'cell junction assembly',
        'cell junction organization',
        'intercellular bridge organization',
        'gap junction-mediated intercellular transport',

        # ECM
        'extracellular matrix assembly',
        'extracellular matrix organization',
        'extracellular matrix constituent secretion',
        'extracellular matrix-cell signaling',
        'cell-matrix recognition',
        'collagen metabolic process',

        # motility
        'fibroblast migration',
        'epithelial cell migration',
        'endothelial cell migration',
        'leukocyte migration',
        'substrate-dependent cell migration',
        'cell chemotaxis to fibroblast growth factor',
        'endothelial cell chemotaxis',
        'fibroblast chemotaxis',
        'leukocyte chemotaxis',
        'epithelial structure maintenance',
        'connective tissue replacement', # this is fibrosis actually

        # channels
        'ion channel activity',

        # docking
        'membrane docking',
        'protein to membrane docking',
        'membrane to membrane docking',

        # kidney specific stuff
        'nephron',
        'outer medulla of kidney',
        'inner medulla of kidney',
        'kidney pyramid',

    },

}

excludes = {
    # the proteins below are not receptors,
    # are excluded from all receptor categories
    'receptor':
        {
            'A6NFA1', 'B2RUY7', 'B4DS77', 'O00170', 'O00468', 'O00548',
            'O00555', 'O14493', 'O14638', 'O14672', 'O14788', 'O14795',
            'O15162', 'O15374', 'O15551', 'O15554', 'O43184', 'O43525',
            'O43813', 'O43866', 'O43914', 'O43921', 'O60291', 'O60359',
            'O75078', 'O75106', 'O75144', 'O75508', 'O75575', 'O75923',
            'O94772', 'O94779', 'O94856', 'O95196', 'O95259', 'O95477',
            'O95727', 'O95897', 'O95967', 'O95994', 'O95998', 'P00734',
            'P00747', 'P01112', 'P01133', 'P01303', 'P01889', 'P01903',
            'P01906', 'P01909', 'P01911', 'P01920', 'P02649', 'P04004',
            'P04439', 'P05026', 'P05156', 'P05187', 'P05231', 'P05362',
            'P05538', 'P07942', 'P08034', 'Q9Y6Y9', 'P09326', 'P09917',
            'P09923', 'P10321', 'P10589', 'P11168', 'P12830', 'P13056',
            'P13569', 'P13591', 'P13598', 'P13637', 'P13765', 'P15514',
            'P15813', 'P16422', 'P16581', 'P17302', 'P17693', 'P17813',
            'P18433', 'P19256', 'P19801', 'P20827', 'P20916', 'P21246',
            'P21589', 'P21926', 'P22001', 'P22460', 'P22736', 'P23276',
            'P23510', 'P23515', 'P24043', 'P25098', 'P25189', 'P26038',
            'P27701', 'P27824', 'P28906', 'P29016', 'P29033', 'P29460',
            'P29972', 'P30301', 'P30511', 'P30533', 'P31431', 'P31749',
            'P31997', 'P32004', 'P32942', 'P32970', 'P32971', 'P35212',
            'P35499', 'P36021', 'P36269', 'P36383', 'P41235', 'P48050',
            'P48509', 'P48552', 'P48651', 'P50591', 'P51787', 'P51828',
            'P52797', 'P52798', 'P52799', 'P52803', 'P52961', 'P54709',
            'P54750', 'P54851', 'P55157', 'P55160', 'P55268', 'P56705',
            'P57057', 'P57087', 'P57739', 'P58401', 'P61160', 'P61769',
            'P62079', 'P62955', 'P78504', 'P78509', 'P78536', 'P78562',
            'P79483', 'P84022', 'P84077', 'P84157', 'P98172', 'Q00994',
            'Q01064', 'Q02108', 'Q02153', 'Q02246', 'Q02413', 'Q02641',
            'Q02846', 'Q03135', 'Q04724', 'Q05940', 'Q06643', 'Q07075',
            'Q07326', 'Q08380', 'Q08AM6', 'Q10588', 'Q12809', 'Q12884',
            'Q12933', 'Q13061', 'Q13114', 'Q13275', 'Q13445', 'Q13520',
            'Q13740', 'Q13822', 'Q13936', 'Q14123', 'Q14126', 'Q14213',
            'Q14242', 'Q14524', 'Q14541', 'Q14563', 'Q14773', 'Q14956',
            'Q14982', 'Q14994', 'Q14995', 'Q15125', 'Q15628', 'Q15758',
            'Q15768', 'Q15842', 'Q15848', 'Q16342', 'Q16363', 'Q16572',
            'Q16625', 'Q16853', 'Q19T08', 'Q29983', 'Q2TAL6', 'Q30201',
            'Q401N2', 'Q496H8', 'Q4VCS5', 'Q4VX76', 'Q4W5P6', 'Q5DX21',
            'Q5T442', 'Q5T4B2', 'Q5TAT6', 'Q5VU97', 'Q5VY80', 'Q5ZPR3',
            'Q6NYC1', 'Q6PIZ9', 'Q6Q4G3', 'Q6RW13', 'Q6UW88', 'Q6UWV6',
            'Q6UXB3', 'Q6UXI9', 'Q7L0J3', 'Q7Z7D3', 'Q86UF1', 'Q86UK0',
            'Q86UR5', 'Q8IU54', 'Q8IU80', 'Q8IUK5', 'Q8IWV2', 'Q8IZV2',
            'Q8IZY2', 'Q8N126', 'Q8N2G4', 'Q8NCM2', 'Q8NCS7', 'Q8NEC5',
            'Q8NFK1', 'Q8NFP4', 'Q8NFY4', 'Q8NFZ3', 'Q8NG11', 'Q8TAZ6',
            'Q8TCY5', 'Q8TD07', 'Q8TDM5', 'Q8WUM9', 'Q8WWA0', 'Q8WWX8',
            'Q8WXS5', 'Q92570', 'Q92673', 'Q92753', 'Q92913', 'Q92954',
            'Q96AP7', 'Q96B86', 'Q96DZ9', 'Q96JB6', 'Q96JK4', 'Q96JQ0',
            'Q96L42', 'Q96PB7', 'Q96QT4', 'Q96S97', 'Q96SN7', 'Q99075',
            'Q99523', 'Q99712', 'Q99784', 'Q99965', 'Q9BQQ7', 'Q9BRK0',
            'Q9BUF7', 'Q9BX67', 'Q9BXJ0', 'Q9BY67', 'Q9BYE2', 'Q9BZM4',
            'Q9GZM7', 'Q9H0R3', 'Q9H221', 'Q9H222', 'Q9H4B8', 'Q9H7V2',
            'Q9NP59', 'Q9NQC3', 'Q9NQS3', 'Q9NR61', 'Q9NR82', 'Q9NRQ2',
            'Q9NSA2', 'Q9NY72', 'Q9NY84', 'Q9NYJ7', 'Q9NYZ4', 'Q9NZ08',
            'Q9NZ53', 'Q9NZQ7', 'Q9NZV8', 'Q9P0K1', 'Q9P0L9', 'Q9P232',
            'Q9P2K9', 'Q9P2U7', 'Q9UBH0', 'Q9UBN1', 'Q9UBX5', 'Q9UEF7',
            'Q9UF02', 'Q9UGM3', 'Q9UHC9', 'Q9UIR0', 'Q9UJZ1', 'Q9UKR5',
            'Q9UKV5', 'Q9UKY0', 'Q9UL54', 'Q9ULB1', 'Q9ULD8', 'Q9ULH0',
            'Q9ULT6', 'Q9UMD9', 'Q9UMF0', 'Q9UNG2', 'Q9UPU3', 'Q9UQ05',
            'Q9UQD0', 'Q9Y215', 'Q9Y219', 'Q9Y286', 'Q9Y2I2', 'Q9Y2J0',
            'Q9Y3R0', 'Q9Y466', 'Q9Y4C0', 'Q9Y566', 'Q9Y5R2', 'Q9Y5Y6',
            'Q9Y5Y9', 'Q9Y624', 'Q9Y698',
        },
    # the proteins below are not ligands,
    # are excluded from all ligand categories
    'ligand':
        {
            'O00220', 'O00300', 'O00468', 'O00587', 'O00592', 'O14594',
            'O14649', 'O14672', 'O15031', 'O15197', 'O15496', 'O43157',
            'O43184', 'O43278', 'O43852', 'O43914', 'O60462', 'O60469',
            'O60486', 'O60494', 'O75051', 'O75077', 'O75330', 'O75508',
            'O75509', 'O75534', 'O75596', 'O94887', 'O95084', 'O95236',
            'O95274', 'O95428', 'O95467', 'O95711', 'O95754', 'P00451',
            'P00488', 'P00533', 'P00734', 'P00740', 'P00742', 'P00748',
            'P00749', 'P00750', 'P00797', 'P00813', 'P00995', 'P01008',
            'P01009', 'P01023', 'P01024', 'P01031', 'P01033', 'P01112',
            'P01589', 'P02452', 'P02458', 'P02461', 'P02462', 'P02654',
            'P02671', 'P02675', 'P02679', 'P02741', 'P02745', 'P02746',
            'P02751', 'P02753', 'P02765', 'P02768', 'P02788', 'P03951',
            'P03956', 'P03973', 'P04003', 'P04070', 'P04196', 'P04278',
            'P04626', 'P04899', 'P05067', 'P05106', 'P05107', 'P05121',
            'P05155', 'P05543', 'P05556', 'P05997', 'P06454', 'P06734',
            'P06756', 'P06858', 'P07093', 'P07225', 'P07237', 'P07288',
            'P07602', 'P07900', 'P08034', 'P08069', 'P08123', 'P08174',
            'P08236', 'P08253', 'P08514', 'P08571', 'P08572', 'P08575',
            'P08582', 'P08603', 'P08670', 'P08709', 'P08861', 'P08865',
            'P09211', 'P09237', 'P09417', 'P09429', 'P0C0L4', 'P0C7T3',
            'P0CG37', 'P0DMV8', 'P0DP23', 'P0DP24', 'P0DP25', 'P10144',
            'P10153', 'P10586', 'P10646', 'P11150', 'P11226', 'P11229',
            'P11362', 'P11912', 'P12107', 'P12109', 'P12110', 'P12821',
            'P12830', 'P13385', 'P13591', 'P13612', 'P13637', 'P13688',
            'P14416', 'P14618', 'P14625', 'P14780', 'P15151', 'P15531',
            'P15907', 'P16035', 'P16070', 'P16109', 'P16473', 'P16520',
            'P16671', 'P17302', 'P17752', 'P17948', 'P18564', 'P19113',
            'P19823', 'P19835', 'P20062', 'P20273', 'P20292', 'P20309',
            'P20849', 'P20908', 'P20916', 'P21462', 'P21709', 'P21802',
            'P21810', 'P21815', 'P21860', 'P21917', 'P21941', 'P21980',
            'P22392', 'P22692', 'P22897', 'P23471', 'P23515', 'P25063',
            'P25090', 'P25940', 'P25942', 'P26012', 'P26441', 'P26842',
            'P27658', 'P28335', 'P28907', 'P28908', 'P29317', 'P29322',
            'P29323', 'P29400', 'P30530', 'P30533', 'P30542', 'P30874',
            'P31025', 'P31785', 'P35212', 'P35354', 'P35475', 'P35555',
            'P35613', 'P35625', 'P36383', 'P36941', 'P39019', 'P39060',
            'P39900', 'P41143', 'P41594', 'P42081', 'P42127', 'P43121',
            'P43405', 'P43489', 'P43490', 'P45452', 'P48039', 'P48651',
            'P49768', 'P49913', 'P50052', 'P50552', 'P51654', 'P52945',
            'P53420', 'P54577', 'P54753', 'P54756', 'P54760', 'P54762',
            'P54764', 'P55058', 'P55789', 'P56159', 'P61626', 'P61769',
            'P62987', 'P63092', 'P78310', 'P78324', 'P78536', 'P80188',
            'P80303', 'P84077', 'P84996', 'P98160', 'P98164', 'Q01469',
            'Q01955', 'Q02094', 'Q02388', 'Q02817', 'Q04721', 'Q05707',
            'Q07326', 'Q07954', 'Q08334', 'Q08722', 'Q08828', 'Q10588',
            'Q12913', 'Q12918', 'Q12933', 'Q13158', 'Q13241', 'Q13255',
            'Q13352', 'Q13361', 'Q13442', 'Q13443', 'Q13444', 'Q13477',
            'Q13936', 'Q14031', 'Q14050', 'Q14055', 'Q14118', 'Q14242',
            'Q14766', 'Q15165', 'Q15223', 'Q15262', 'Q15303', 'Q15375',
            'Q15762', 'Q16613', 'Q2MV58', 'Q2VPA4', 'Q4VX76', 'Q5JWF2',
            'Q5SR53', 'Q5T442', 'Q5T5A4', 'Q6NW40', 'Q6UWW8', 'Q6UWX4',
            'Q6V0I7', 'Q7Z6A9', 'Q86UR5', 'Q8IWL1', 'Q8IZJ3', 'Q8IZL2',
            'Q8N2X6', 'Q8N474', 'Q8NFK1', 'Q8NFT8', 'Q8NGH5', 'Q8NGH8',
            'Q8NHJ6', 'Q8NHP8', 'Q8TAX7', 'Q8WTV0', 'Q8WWY8', 'Q8WZ79',
            'Q92819', 'Q92854', 'Q92896', 'Q92956', 'Q96A49', 'Q96CG8',
            'Q96DA0', 'Q96FE5', 'Q99259', 'Q99466', 'Q99965', 'Q9BQ66',
            'Q9BS26', 'Q9BX66', 'Q9BZR6', 'Q9BZW8', 'Q9BZZ2', 'Q9C0C4',
            'Q9H2A7', 'Q9H2E6', 'Q9H3S1', 'Q9H9H4', 'Q9HCM2', 'Q9NQC3',
            'Q9NR96', 'Q9NRV9', 'Q9NSG2', 'Q9NTN9', 'Q9NUP9', 'Q9NV23',
            'Q9NWZ3', 'Q9NX52', 'Q9NZC2', 'Q9NZR2', 'Q9UBY5', 'Q9UHG3',
            'Q9UIW2', 'Q9UJU6', 'Q9UKQ2', 'Q9ULP9', 'Q9UM47', 'Q9UQ26',
            'Q9XRX5', 'Q9Y215', 'Q9Y2I2', 'Q9Y4D7', 'Q9Y566', 'Q9Y5U5',
            'Q9Y624', 'Q9Y625', 'Q9Y6N7', 'P11362', 'P14778', 'P36897',
        },
    # the proteins below are not cell surface ligand proteins,
    # are excluded from all cell surface ligand categories
    'cell_surface_ligand': {
            'P00533', 'P0DPD6',
        },
    # the proteins below are not adhesion proteins,
    # are excluded from all adhesion categories
    'adhesion':
        {
            'O14976', 'O15231', 'O43294', 'O60711', 'O60759', 'P07332',
            'P29965', 'P35221', 'P35222', 'P50552', 'P56945', 'Q02297',
            'Q03001', 'Q13153', 'Q13895', 'Q14451', 'Q14511', 'Q14943',
            'Q14952', 'Q14953', 'Q155Q3', 'Q15654', 'Q5JRA6', 'Q5T4B2',
            'Q7L5Y9', 'Q7Z4I7', 'Q86TP1', 'Q8IVT2', 'Q8IZW8', 'Q8N264',
            'Q8N743', 'Q8NHK3', 'Q8WX93', 'Q92502', 'Q96AC1', 'Q96IF1',
            'Q96QB1', 'Q99689', 'Q9H792', 'Q9HBI0', 'Q9HBI1', 'Q9NQ75',
            'Q9UBT7', 'Q9UGI8', 'Q9UGP4', 'Q9UI47', 'Q9UQB3',
        },
    # the proteins below are not cell surface enzymes,
    # are excluded from all cell surface enzyme categories
    'cell_surface_enzyme':
        {
            'P52961', 'O00391', 'Q8WUD6', 'P11117', 'Q15125', 'Q13444',
            'P08842', 'Q9UK23', 'Q9NPH5', 'P19021', 'P09958', 'P40126',
            'Q96JJ7', 'P04843', 'Q8TCJ2', 'Q9UKF2',
        },
}


"""
Higher level classes of intercellular communication roles.
"""
annot_combined_classes = (

    ### locational (and topological) categories ###

    # transmembrane
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        source = 'composite',
        scope = 'generic',
        resource = '~transmembrane',
    ),
    af.AnnotDef(
        name = 'transmembrane',
        parent = 'transmembrane',
        resource = 'UniProt_location',
        aspect = 'locational',
        scope = 'generic',
        args = {
            'features': {
                'Multi-pass membrane protein',
                'Single-pass membrane protein',
                'Single-pass type I membrane protein',
                'Single-pass type II membrane protein',
                'Single-pass type III membrane protein',
                'Single-pass type IV membrane protein',
            },
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        parent = 'transmembrane',
        scope = 'generic',
        resource = 'UniProt_topology',
        aspect = 'locational',
        args = {
            'topology': 'Transmembrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        parent = 'transmembrane',
        scope = 'generic',
        resource = 'UniProt_keyword',
        aspect = 'locational',
        args = {
            'keyword': {
                'Transmembrane',
                'Transmembrane beta strand',
                'Transmembrane helix',
            },
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_predicted',
        parent = 'transmembrane',
        scope = 'generic',
        source = 'composite',
        aspect = 'locational',
        resource = af.AnnotOp(
            annots = '~transmembrane_predicted',
            op = common.at_least_in(2),
        ),
    ),
    af.AnnotDef(
        name = 'transmembrane',
        parent = 'transmembrane_predicted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Phobius',
        args = {
            'tm_helices': bool,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_phobius',
        parent = 'transmembrane_predicted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Almen2009',
        args = {
            'phobius_transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_sosui',
        parent = 'transmembrane_predicted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Almen2009',
        args = {
            'sosui_transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane_tmhmm',
        parent = 'transmembrane_predicted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Almen2009',
        args = {
            'tmhmm_transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'transmembrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        resource = 'CellPhoneDB',
        args = {
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'OPM',
        args = {
            'transmembrane': True,
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'TopDB',
        args = {
            'topology': 'Membrane',
        },
    ),
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'LOCATE',
        args = {
            'cls': {
                'typeI',
                'typeII',
                'mtmp',
            },
        },
    ),  # about 60 proteins above the ones in UniProt classified as
        # transmembrane, mostly based on prediction methods
        # for this reason we don't use it, however it might be that
        # many more proteins have transmembrane isoforms apart from
        # the ones in UniProt
    af.AnnotDef(
        name = 'transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'tmh': bool,
        },
    ),  # same as for LOCATE: overall 141 additional TM proteins
        # apart from the ones in UniProt
    af.AnnotDef(
        name = 'lhfpl',
        parent = 'transmembrane',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'LHFPL tetraspan proteins',
        },
    ),

    # peripheral
    af.AnnotDef(
        name = 'peripheral',
        source = 'composite',
        scope = 'generic',
        aspect = 'locational',
        resource = '~peripheral',
    ),
    af.AnnotDef(
        name = 'peripheral',
        parent = 'peripheral',
        scope = 'generic',
        resource_name = 'UniProt_location',
        aspect = 'locational',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'peripheral_lipid_anchor',
                    resource = 'UniProt_location',
                    args = {
                        'features': {
                            'Peripheral membrane protein',
                            'Lipid-anchor',
                        },
                    },
                ),
                af.AnnotDef(
                    name = 'peripheral_gpi_anchor',
                    resource = 'UniProt_location',
                    args = {
                        'location': {
                            'GPI-anchor',
                            'GPI-like-anchor',
                        },
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'peripheral',
        parent = 'peripheral',
        scope = 'generic',
        resource = 'UniProt_topology',
        aspect = 'locational',
        args = {
            'topology': 'Intramembrane',
        },
    ),

    # plasma membrane
    af.AnnotDef(
        name = 'plasma_membrane',
        source = 'composite',
        aspect = 'locational',
        scope = 'generic',
        resource = '~plasma_membrane',
    ),
    af.AnnotDef(
        name = 'plasma_membrane',
        parent = 'plasma_membrane',
        scope = 'generic',
        resource = 'UniProt_location',
        aspect = 'locational',
        args = {
            'location': {
                'Cell membrane',
                'Basal cell membrane',
                'Basolateral cell membrane',
                'Lateral cell membrane',
                'Apicolateral cell membrane',
                'Apical cell membrane',
            },
        },
    ),
    af.AnnotDef(
        name = 'plasma_membrane',
        parent = 'plasma_membrane',
        scope = 'generic',
        resource = 'Cellinker',
        aspect = 'locational',
        args = {
            'location': 'Membrane'
        },
    ),

    # plasma membrane regions
    # from UniProt_location
    af.AnnotDef(
        name = 'basolateral_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Basolateral cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'basal_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Basal cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'apical_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Apical cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'apicolateral_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Apicolateral cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'lateral_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Lateral cell membrane',
        },
    ),
    # from Ramilowski_location
    af.AnnotDef(
        name = 'plasma_membrane',
        resource = '~plasma_membrane~Ramilowski_location',
        resource_name = 'Ramilowski_location',
        scope = 'generic',
        aspect = 'locational',
    ),
    af.AnnotDef(
        name = 'basolateral_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'location': 'basolateral cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'basal_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'location': 'basal cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'apical_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'location': 'apical cell membrane',
        },
    ),
    af.AnnotDef(
        name = 'lateral_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'location': 'lateral cell membrane',
        },
    ),
    # from LOCATE
    af.AnnotDef(
        name = 'plasma_membrane',
        resource = '~plasma_membrane~LOCATE',
        resource_name = 'LOCATE',
        scope = 'generic',
        aspect = 'locational',
    ),
    af.AnnotDef(
        name = 'basolateral_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'LOCATE',
        args = {
            'location': 'basolateral plasma membrane',
        },
    ),
    af.AnnotDef(
        name = 'apical_cell_membrane',
        parent = 'plasma_membrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'LOCATE',
        args = {
            'location': 'apical plasma membrane',
        },
    ),

    # plasma membrane transmembrane
    af.AnnotDef(
        name = 'plasma_membrane_transmembrane',
        source = 'composite',
        scope = 'generic',
        aspect = 'locational',
        resource = af.AnnotOp(
            annots = (
                af.AnnotOp(
                    annots = (
                        'transmembrane',
                        'plasma_membrane',
                    ),
                    op = set.intersection,
                ),
                '~plasma_membrane_transmembrane',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'plasma_membrane_transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'Membranome',
        args = {
            'membrane': 'Plasma membrane',
            'side': 'extracellular side',
        },
        exclude = {
            'O14798', 'O75326', 'P04216', 'Q6H3X3', 'P55259', 'P22748',
        },
    ),  # with a few exception these are transmembrane proteins
        # of the plasma membrane
    af.AnnotDef(
        name = 'plasma_membrane_transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'CSPA',
        args = {
            'high_confidence': bool,
            'tm': bool,
        },
    ),
    af.AnnotDef(
        name = 'plasma_membrane_transmembrane',
        aspect = 'locational',
        scope = 'generic',
        resource = 'HPMR',
        args = {
            'role': 'Receptor',
        },
        exclude = {
            'P56159', 'Q14982', 'P35052', 'O14798', 'Q96QV1', 'P26992',
            'Q9Y5V3', 'O00451', 'Q12860', 'O60609', 'P14207', 'O43813',
            'P15328', 'O75015', 'Q9BZR6',
        },
    ),
    af.AnnotDef(
        name = 'ifn_induced',
        parent = 'plasma_membrane_transmembrane',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon induced transmembrane proteins',
        },
    ),
    af.AnnotDef(
        name = 'lhfpl_plasma_membrane',
        parent = 'plasma_membrane_transmembrane',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'LHFPL tetraspan proteins',
        }
    ),

    #### what to do with this???
    # plasma membrane regulator
    af.AnnotDef(
        name = 'plasma_membrane_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'EMP-PMP22-LIM',
        },
    ),

    # plasma membrane peripheral
    af.AnnotDef(
        name = 'plasma_membrane_peripheral',
        source = 'composite',
        scope = 'generic',
        aspect = 'locational',
        resource = af.AnnotOp(
            annots = (
                af.AnnotOp(
                    annots = (
                        'peripheral',
                        'plasma_membrane',
                    ),
                    op = set.intersection,
                ),
                '~plasma_membrane_peripheral',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'plasma_membrane_peripheral',
        aspect = 'locational',
        scope = 'generic',
        resource = 'CSPA',
        args = {
            'high_confidence': bool,
            'gpi': bool,
        },
    ),

    # secreted
    af.AnnotDef(
        name = 'secreted',
        source = 'composite',
        scope = 'generic',
        aspect = 'locational',
        resource = '~secreted',
    ),
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        resource = 'UniProt_keyword',
        aspect = 'locational',
        args = {
            'keyword': 'Secreted',
        },
    ),
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        resource = 'UniProt_location',
        aspect = 'locational',
        args = {
            'location': 'Secreted',
        },
    ),
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        resource = 'HPA_secretome',
        aspect = 'locational',
        args = {
            'secreted': bool,
        },
    ),  # looks all right
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'connectomeDB2020',
        args = {
            'location': 'secreted',
        },
    ),
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        resource = 'Baccin2019',
        aspect = 'locational',
        args = {
            'location': {'secreted', 'both', 'ecm'},
        },
    ),
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        resource = 'MatrixDB',
        aspect = 'locational',
        args = {
            'mainclass': 'secreted',
        },
        enabled = False,
    ),  # some potentially wrong elements, proteins annotated by
        # UniProt as intracellular
        # manual revision is time consuming (3000 proteins), disabled
        # for the time being
    af.AnnotDef(
        name = 'secreted',
        aspect = 'locational',
        scope = 'generic',
        resource = 'LOCATE',
        args = {
            'cls': 'secretome',
        },
        enabled = False,
    ),  # unusable, too many intracellular proteins
        # secreted
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Matrisome',
        args = {
            'subclass': 'Secreted Factors',
        },
        exclude = {'P51610'}
    ),
    af.AnnotDef(
        name = 'secreted',
        parent = 'secreted',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Cellinker',
        args = {
            'location': 'Secreted',
        },
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'bpi_fold_containing',
        parent = 'secreted',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'BPI fold containing',
        },
    ),
    af.AnnotDef(
        name = 'histatin',
        parent = 'secreted',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Histatins and statherin',
        },
    ),  # secreted into saliva
    af.AnnotDef(
        name = 'proline_rich',
        parent = 'secreted',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Proline rich proteins',
        },
    ),  # secreted into saliva, enamel protective and anti-microbial

    # cell_surface =
    # plasma_membrane_peripheral + plasma_membrane_transmembrane
    af.AnnotDef(
        name = 'cell_surface',
        parent = 'cell_surface',
        scope = 'generic',
        aspect = 'locational',
        source = 'composite',
        resource = af.AnnotOp(
            annots = (
                'plasma_membrane_transmembrane',
                'plasma_membrane_peripheral',
                '~cell_surface',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'cell_surface',
        parent = 'cell_surface',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Surfaceome',
        exclude = {
            'Q7L1I2', 'Q9ULQ1', 'Q05940', 'Q9BZC7', 'Q8NBW4', 'P54219',
            'Q9P2U8', 'Q8IY34', 'Q8TED4', 'Q9UN42', 'Q9P2U7', 'Q8NCC5',
            'Q9H598', 'Q8NHS3', 'Q9NRX5', 'Q9H1V8', 'Q496J9', 'Q6J4K2',
            'Q96T83', 'Q9NP78', 'A6NFC5', 'Q8TBB6', 'O00400', 'Q8WWZ7',
            'Q71RS6', 'Q9GZU1', 'O95528', 'Q8NDX2', 'O43826', 'O94778',
            'Q9HD20', 'Q9UGQ3', 'Q14108',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface',
        parent = 'cell_surface',
        scope = 'generic',
        aspect = 'locational',
        resource = 'connectomeDB2020',
        args = {
            'location': 'plasma membrane',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface',
        parent = 'cell_surface',
        scope = 'generic',
        resource = 'Baccin2019',
        aspect = 'locational',
        args = {
            'location': {'membrane', 'both'},
        },
    ),
    af.AnnotDef(
        name = 'cell_surface',
        resource = 'Ramilowski_location',
        scope = 'generic',
        aspect = 'locational',
        args = {
            'location': 'cell surface',
        },
        enabled = False,
    ),  # mostly intracellular, disabled
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'glypican',
        parent = 'cell_surface',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glypicans',
        },
    ),
    af.AnnotDef(
        name = 'glypican',
        parent = 'cell_surface',
        aspect = 'locational',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Glypican',
        },
    ),
    af.AnnotDef(
        name = 'syndecan',
        parent = 'cell_surface',
        aspect = 'locational',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Syndecan',
        },
    ),
    af.AnnotDef(
        name = 'immunoglobulin_heavy',
        parent = 'cell_surface',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin heavy locus at 14q32.33',
        },
    ),  # immunoglobulin heavy chain
    af.AnnotDef(
        name = 'immunoglobulin_kappa',
        parent = 'cell_surface',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin kappa locus at 2p11.2',
        },
    ),  # immunoglobulin V region
    af.AnnotDef(
        name = 'immunoglobulin_lambda',
        parent = 'cell_surface',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin lambda locus at 22q11.2',
        },
    ),  # immunoglobulin V region

    # extracellular =
    # cell_surface + secreted + ecm
    af.AnnotDef(
        name = 'extracellular',
        parent = 'extracellular',
        scope = 'generic',
        aspect = 'locational',
        source = 'composite',
        resource = af.AnnotOp(
            annots = (
                'secreted',
                'cell_surface',
                'ecm',
                '~extracellular',
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'extracellular',
        parent = 'extracellular',
        scope = 'generic',
        aspect = 'locational',
        resource = 'Ramilowski_location',
        args = {
            'location': {
                'secreted',
            },
        },
        enabled = False,
    ),  # many membrane bound, transmembrane or intracellular
        # at least according to UniProt
        # disabled because of this
    af.AnnotDef(
        name = 'extracellular',
        parent = 'extracellular',
        scope = 'generic',
        aspect = 'locational',
        resource = 'HPMR',
    ),  # these are membrane bound or secreted proteins
        # so we can add them only to the extracellular
    af.AnnotDef(
        name = 'extracellular',
        resource = 'DGIdb',
        aspect = 'locational',
        args = {
            'category': 'EXTERNAL SIDE OF PLASMA MEMBRANE',
        },
        exclude = {
            'O75534'
        },
    ),  # the `CELL SURFACE` category is completely unusable, contains
        # proteins from any location randomly
        # this on, the `EXT. SIDE OF PM.` is better, but contains many
        # secreted proteins, ligands, and some potentially intracellular ones
    af.AnnotDef(
        name = 'extracellular',
        aspect = 'locational',
        scope = 'generic',
        resource = 'ComPPI',
        args = {
            'location': 'extracellular',
        },
        enabled = False,
    ),  # unusable, unrealistically huge and
        # contains half of the intracellular proteome
    af.AnnotDef(
        name = 'extracellular',
        scope = 'generic',
        aspect = 'locational',
        resource = 'LOCATE',
        args = {
            'location': {
                'extracellular',
                'extracellular region',
            },
        },
    ),  # unusable, too many intracellular proteins
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'iglon5',
        parent = 'extracellular',
        aspect = 'locational',
        resource = {'A6NGN9'},
    ),  # function is not clear for me
    af.AnnotDef(
        name = 'immunoglobulin_heavy',
        parent = 'extracellular',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin heavy locus at 14q32.33',
        },
    ),  # immunoglobulin heavy chain
    af.AnnotDef(
        name = 'immunoglobulin_kappa',
        parent = 'extracellular',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin kappa locus at 2p11.2',
        },
    ),  # immunoglobulin V region
    af.AnnotDef(
        name = 'immunoglobulin_lambda',
        parent = 'extracellular',
        aspect = 'locational',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin lambda locus at 22q11.2',
        },
    ),  # immunoglobulin V region

    # intracellular
    af.AnnotDef(
        name = 'intracellular',
        aspect = 'locational',
        resource = '~intracellular',
        source = 'composite',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'intracellular',
        aspect = 'locational',
        resource_name = 'LOCATE',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'intracellular',
                    resource = 'LOCATE',
                    args = {
                        'location': {
                            'centrosome',
                            'cytoplasm',
                            'endosomes',
                            'lysosomes',
                            'nucleus',
                            'plasma membrane',
                            'cytoplasmic membrane-bound vesicle',
                            'cytoplasmic vesicles',
                            'cytoskeleton',
                            'early endosomes',
                            'endoplasmic reticulum',
                            'golgi apparatus',
                            'er-golgi intermediate compartment',
                            'ergic',
                            'golgi cis cisterna',
                            'golgi medial cisterna',
                            'golgi trans cisterna',
                            'golgi trans face',
                            'inner mitochondrial membrane',
                            'late endosomes',
                            'lipid particles',
                            'medial-golgi',
                            'melanosome',
                            'microtubule',
                            'microtubule organizing center ',
                            'mitochondria',
                            'mitochondrial inner membrane',
                            'mitochondrial outer membrane',
                            'mitochondrion',
                            'nuclear envelope',
                            'nucleolus',
                            'nuclear speck',
                            'outer mitochondrial membrane',
                            'peroxisome',
                            'peroxisomes',
                            'sarcolemma',
                            'transport vesicle',
                        },
                    },
                ),
                af.AnnotDef(
                    name = 'cytoplasmic',
                    resource = 'LOCATE',
                    args = {
                        'cls': 'cytoplasmic',
                    },
                ),
            ),
            op = set.union,
        ),
    ),
    af.AnnotDef(
        name = 'intracellular',
        aspect = 'locational',
        scope = 'generic',
        resource = 'ComPPI',
        args = {
            'location': {
                'cytosol',
                'nucleus',
                'mitochondrion',
            },
        },
    ),
    af.AnnotDef(
        name = 'intracellular',
        aspect = 'locational',
        scope = 'generic',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'intracellular',
        },
    ),
    af.AnnotDef(
        name = 'intracellular',
        aspect = 'locational',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': {
                'Autophagosome',
                'Autophagosome lumen',
                'Autophagosome membrane',
                'Centriolar satellite',
                'Centriole',
                'Centromere',
                'Centrosome',
                'Chromaffin granule membrane',
                'Chromosome',
                'Cilium basal body',
                'Cis-Golgi network',
                'Cis-Golgi network membrane',
                'Cytoplasm',
                'Cytoplasmic granule',
                'Cytoplasmic granule lumen',
                'Cytoplasmic granule membrane',
                'Cytoplasmic vesicle',
                'Cytoplasmic vesicle membrane',
                'Cytoskeleton',
                'Cytosol',
                'Early endosome',
                'Early endosome membrane',
                'Endomembrane system',
                'Endoplasmic reticulum',
                'Endoplasmic reticulum lumen',
                'Endoplasmic reticulum membrane',
                'Endoplasmic reticulum-Golgi intermediate compartment', (
                    'Endoplasmic reticulum-Golgi intermediate '
                    'compartment membrane'
                ),
                'Endosome',
                'Endosome lumen',
                'Endosome membrane',
                'Golgi apparatus',
                'Golgi apparatus lumen',
                'Golgi apparatus membrane',
                'Golgi stack',
                'Golgi stack membrane',
                'Late endosome',
                'Late endosome membrane',
                'Lysosome',
                'Lysosome lumen',
                'Lysosome membrane',
                'Melanosome',
                'Melanosome membrane',
                'Microsome',
                'Microsome membrane',
                'Microtubule organizing center',
                'Mitochondrion',
                'Mitochondrion inner membrane',
                'Mitochondrion intermembrane space',
                'Mitochondrion matrix',
                'Mitochondrion membrane',
                'Mitochondrion nucleoid',
                'Mitochondrion outer membrane',
                'Multivesicular body',
                'Nuclear pore complex',
                'Nucleolus',
                'Nucleoplasm',
                'Nucleus',
                'Nucleus envelope',
                'Nucleus inner membrane',
                'Nucleus lamina',
                'Nucleus matrix',
                'Nucleus membrane',
                'Nucleus outer membrane',
                'Nucleus speckle',
                'P-body',
                'PML body',
                'Perikaryon',
                'Perinuclear region',
                'Perinuclear theca',
                'Peroxisome',
                'Peroxisome matrix',
                'Peroxisome membrane',
                'Postsynaptic Golgi apparatus',
                'Preautophagosomal structure',
                'Preautophagosomal structure membrane',
                'Recycling endosome',
                'Recycling endosome membrane',
                'Rough endoplasmic reticulum',
                'Rough endoplasmic reticulum lumen',
                'Rough endoplasmic reticulum membrane',
                'Sarcomere',
                'Sarcoplasmic reticulum',
                'Sarcoplasmic reticulum lumen',
                'Sarcoplasmic reticulum membrane',
                'Smooth endoplasmic reticulum membrane',
                'Telomere',
                'Trans-Golgi network',
            },
        },
    ),

    ### functional classes ###

    # receptor
    af.AnnotDef(
        name = 'receptor',
        resource = '~receptor',
        scope = 'generic',
        source = 'composite',
        receiver = True,
        transmitter = False,
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'talklr',
        args = {
            'role': 'receptor',
            'putative': False,
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'Cellinker',
        args = {
            'role': 'receptor',
            'type': {
                'Cytokine-cytokine receptor interaction',
                'Secreted protein to receptor interaction',
            },
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'scConnect',
        args = {
            'role': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'connectomeDB2020',
        args = {
            'role': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'CellCall',
        args = {
            'role': 'receptor',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'iTALK',
        args = {
            'mainclass': 'receptor',
        },
        scope = 'generic',
        exclude = {
            'Q13114', 'P23510', 'Q9NZV8', 'P48050', 'Q9NSA2', 'P22001',
            'Q9NR82', 'O15554', 'P51787', 'O43525', 'Q9UKR5', 'Q99712',
            'P78562', 'Q9UQD0', 'Q9Y5Y9', 'Q9ULH0', 'Q03135', 'P52961',
            'Q14524', 'P35499', 'P27701', 'Q00994', 'Q9NZ08', 'Q9UMD9',
            'Q05940', 'Q16572', 'Q9P2U7', 'Q8NCS7', 'O15374', 'P36021',
            'Q8WWX8', 'P57057', 'Q9Y5R2', 'Q96SN7', 'O95994', 'Q14123',
            'Q01064', 'P54750', 'Q9NQS3', 'P29033', 'Q5VU97', 'Q96JQ0',
            'O43914', 'Q96AP7', 'Q16625', 'P13637', 'P09917', 'P57739',
            'P51828', 'Q13936', 'P29016', 'Q13520', 'O60291', 'P30301',
            'Q8NEC5', 'O95477', 'O14493', 'P11279', 'P09326', 'P11168',
            'P16581', 'P18433', 'P20916', 'P21926', 'P29972', 'P30511',
            'P31431', 'P48509', 'P58401', 'Q02846', 'Q12884', 'Q30201',
            'Q7Z7D3', 'Q9NP59', 'Q9NRQ2', 'Q9ULB1', 'Q9Y4C0', 'Q9Y624',
            'Q9Y6Y9', 'Q12933', 'Q5TAT6', 'P34741', 'P78508', 'Q9HCN3',
            'Q8N158', 'Q5TAH2', 'Q9UHC6', 'Q9ULT6', 'Q5HYA8', 'P40145',
            'P21589', 'Q15746',
        },
    ),  # locations are correct, includes also adhesion receptors
    af.AnnotDef(
        name = 'growth_factor',
        parent = 'receptor',
        resource = 'iTALK',
        args = {
            'mainclass': 'receptor',
            'subclass': 'growth factor',
        },
    ),  # looks good
    af.AnnotDef(
        name = 'cytokine',
        parent = 'receptor',
        resource = 'iTALK',
        args = {
            'mainclass': 'receptor',
            'subclass': 'cytokine',
        },
    ),  # looks good
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'CellCellInteractions',
        args = {
            'mainclass': 'Receptor',
        },
        exclude = {
            'Q8IZV2', 'O95477', 'A6NFA1', 'B2RUY7', 'B4DS77', 'O00170',
            'O00468', 'O00548', 'O00555', 'O14493', 'O14638', 'O14672',
            'O14788', 'O14795', 'O15551', 'O43813', 'O43866', 'O43921',
            'O60359', 'O75078', 'O75106', 'O75144', 'O75575', 'O94772',
            'O94779', 'O94856', 'O95196', 'O95259', 'O95727', 'O95897',
            'O95967', 'P00734', 'P01303', 'P01889', 'P01903', 'P01906',
            'P01909', 'P01911', 'P01920', 'P02649', 'P04004', 'P04439',
            'P05026', 'P05156', 'P05187', 'P05231', 'P05362', 'P05538',
            'Q9P2K9', 'P09326', 'P09923', 'P10321', 'P10589', 'P11168',
            'P12830', 'P13056', 'P13598', 'P13765', 'P15813', 'P16581',
            'P17693', 'P18433', 'P19256', 'P19801', 'P20827', 'P20916',
            'P21246', 'P21926', 'P22460', 'P22736', 'P23515', 'P25098',
            'P29460', 'P29972', 'P30511', 'P30533', 'P31431', 'P31997',
            'P32004', 'P32942', 'P32970', 'P36269', 'P41235', 'P48509',
            'P52797', 'P52798', 'P52799', 'P52803', 'P54709', 'P54851',
            'P55160', 'P56705', 'P58401', 'P62079', 'P62955', 'P78504',
            'P78509', 'P78536', 'P79483', 'P84157', 'P98172', 'Q02108',
            'Q02153', 'Q02246', 'Q02641', 'Q02846', 'Q08380', 'Q08AM6',
            'Q10588', 'Q12809', 'Q12884', 'Q13061', 'Q13275', 'Q13445',
            'Q13740', 'Q13822', 'Q14213', 'Q14242', 'Q14541', 'Q14563',
            'Q14773', 'Q14956', 'Q14982', 'Q14994', 'Q14995', 'Q15125',
            'Q15758', 'Q15768', 'Q15842', 'Q15848', 'Q16853', 'Q19T08',
            'Q29983', 'Q2TAL6', 'Q30201', 'Q401N2', 'Q496H8', 'Q4VCS5',
            'Q4W5P6', 'Q5DX21', 'Q5T4B2', 'Q5VY80', 'Q5ZPR3', 'Q6PIZ9',
            'Q6Q4G3', 'Q6RW13', 'Q6UWV6', 'Q6UXB3', 'Q6UXI9', 'Q7L0J3',
            'Q7Z7D3', 'Q86UF1', 'Q86UK0', 'Q8IU54', 'Q8IUK5', 'Q8IWV2',
            'Q8IZY2', 'Q8N2G4', 'Q8NCM2', 'Q8NFP4', 'Q8NFY4', 'Q8NFZ3',
            'Q8NG11', 'Q8TAZ6', 'Q8TCY5', 'Q8TD07', 'Q8TDM5', 'Q8WUM9',
            'Q8WWA0', 'Q92570', 'Q92753', 'Q92913', 'Q92954', 'Q96B86',
            'Q96DZ9', 'Q96JB6', 'Q96JK4', 'Q96L42', 'Q96PB7', 'Q96S97',
            'Q99075', 'Q99784', 'Q9BQQ7', 'Q9BRK0', 'Q9BUF7', 'Q9BY67',
            'Q9BYE2', 'Q9BZM4', 'Q9GZM7', 'Q9H221', 'Q9H222', 'Q9H4B8',
            'Q9H7V2', 'Q9NP59', 'Q9NR61', 'Q9NRQ2', 'Q9NY72', 'Q9NY84',
            'Q9NYJ7', 'Q9NYZ4', 'Q9P0K1', 'Q9P0L9', 'Q9P232', 'Q9P2K9',
            'Q9UBN1', 'Q9UBX5', 'Q9UEF7', 'Q9UF02', 'Q9UGM3', 'Q9UHC9',
            'Q9UIR0', 'Q9UJZ1', 'Q9UL54', 'Q9ULB1', 'Q9ULD8', 'Q9UMF0',
            'Q9UQ05', 'Q9Y215', 'Q9Y219', 'Q9Y286', 'Q9Y2I2', 'Q9Y3R0',
            'Q9Y466', 'Q9Y4C0', 'Q9Y624', 'Q9Y698', 'Q9Y6Y9', 'Q8WXS5',
            'Q9BX67', 'Q12933', 'Q96QT4', 'P20645',
        },
    ),  # includes also secreted proteins, such as ligands, and intracellular
        # proteins, these we try to exclude;
        # furthermore nuclear receptors or receptors in the endosome
        # membrane, but I think those are OK
    af.AnnotDef(
        name = 'receptor',
        resource = 'EMBRACE',
        args = {
            'mainclass': 'receptor',
        },
        scope = 'generic',
        exclude = {
            'Q30201', 'Q05940', 'Q96SN7', 'P29972', 'Q9NP59', 'O95477',
            'Q9Y624', 'P29033', 'Q5VU97', 'Q9Y6Y9', 'P52961', 'Q9P2U7',
            'P20916', 'Q96JQ0', 'O43914', 'Q96AP7', 'O14493', 'P22001',
            'P09326', 'Q14524', 'Q99712', 'Q9NRQ2', 'Q16625', 'P27701',
            'P13637', 'Q9NZV8', 'P78562', 'P51828', 'P18433', 'Q13936',
            'Q9NSA2', 'P21926', 'Q12933', 'Q12884', 'P31431', 'Q9UQD0',
            'Q03135', 'Q9ULB1', 'Q6NYC1', 'Q8WWX8', 'Q8NEC5', 'P48509',
            'P57057',
        },
    ),  # good contents, includes a few secreted and
        # intracellular receptors
    af.AnnotDef(
        name = 'gpcr',
        parent = 'receptor',
        resource = 'GPCRdb',
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = '~receptor~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'CellPhoneDB',
        args = {
            'receptor': True,
            'transmembrane': True,
        },
        scope = 'generic',
        exclude = {'P14735', 'P52799', 'P05231', 'Q15768', 'Q8IZI9'},
    ),  # includes some cell-cell adhesion
    af.AnnotDef(
        name = 'receptor',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'receptors',
        },
        scope = 'generic',
        exclude = {
            'P11021', 'P19835', 'O75144', 'P35225', 'P56705', 'P12643',
            'Q8IX30', 'Q14512', 'P04628', 'Q76LX8', 'P02647', 'O00755',
            'P04085', 'P10600', 'P14735',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'HPMR',
        args = {
            'role': 'Receptor',
        },
        scope = 'generic',
        # `exclude` defined in the `excludes` dict
    ),  # includes adhesion molecules e.g. intergrins
    af.AnnotDef(
        name = 'receptor',
        resource = 'ICELLNET',
        args = {
            'role': 'receptor',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'CellChatDB',
        args = {
            'role': 'receptor',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'CellTalkDB',
        args = {
            'role': 'receptor',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'Surfaceome',
        args = {
            'mainclass': 'Receptors',
        },
        scope = 'generic',
        exclude = {'Q14108'}
    ),  # good as it is
    af.AnnotDef(
        name = 'receptor',
        resource = 'Ramilowski2015',
        args = {
            'mainclass': 'receptor',
        },
        scope = 'generic',
        exclude = {
            'Q9NZ08', 'Q30201', 'Q9Y5R2', 'Q96SN7', 'O95994', 'P29972',
            'Q9NP59', 'O95477', 'Q9Y624', 'Q9NQS3', 'P52961', 'P20916',
            'P16581', 'Q96JQ0', 'O43914', 'O14493', 'P22001', 'Q14524',
            'Q99712', 'Q9NRQ2', 'P27701', 'Q00994', 'P13637', 'Q9NZV8',
            'P09917', 'Q9Y4C0', 'P57739', 'P78562', 'Q13936', 'P29016',
            'Q9NSA2', 'P21926', 'O60291', 'Q12933', 'P36021', 'Q12884',
            'P31431', 'P30511', 'Q9UQD0', 'Q03135', 'Q02846', 'P48050',
            'P58401', 'Q9ULB1', 'Q7Z7D3', 'Q9UKR5', 'Q08828', 'Q8NEC5',
            'P48509', 'P57057', 'Q5TAT6', 'Q92542', 'Q02094', 'Q9Y6N8',
            'O14939', 'Q9ULT6', 'O75923', 'P21589', 'P13726', 'P13569',
        },
    ),  # many non-receptors excluded
    af.AnnotDef(
        name = 'receptor',
        resource = 'Kirouac2010',
        args = {
            'mainclass': 'receptor',
        },
        scope = 'generic',
    ),  # good as it is
    af.AnnotDef(
        name = 'receptor',
        resource = 'Guide2Pharma',
        args = {
            'mainclass': 'receptor',
        },
        scope = 'generic',
        exclude = {
            'Q9H1R3', 'P54750', 'Q14123', 'P09917', 'Q02846', 'Q01064',
        },
    ),  # these are really true receptors
    af.AnnotDef(
        name = 'gpcr',
        parent = 'receptor',
        resource = 'DGIdb',
        args = {
            'category': 'G PROTEIN COUPLED RECEPTOR',
        },
        exclude = {
            'O43813', 'Q9P2Y4', 'Q6RW13', 'P16112', 'Q9BRX2', 'O75575',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        resource = 'LRdb',
        args = {
            'role': 'receptor',
            'references': bool,
        },
        scope = 'generic',
        exclude = {
            'O14493', 'O15374', 'O15554', 'O43525', 'O43914', 'O60291',
            'O95477', 'O95994', 'O95998', 'P00747', 'Q9UKV5', 'P09326',
            'P09917', 'P11168', 'P13637', 'P15813', 'P16581', 'P18433',
            'P20916', 'P21926', 'P22001', 'P26038', 'P27701', 'P27824',
            'P29016', 'P29033', 'P29972', 'P30301', 'P30511', 'P31431',
            'P35499', 'P36021', 'P48050', 'P48509', 'P51787', 'P51828',
            'P52797', 'P52961', 'P54750', 'P55157', 'P57057', 'P57739',
            'P58401', 'P78536', 'P78562', 'P84022', 'Q00994', 'Q01064',
            'Q02413', 'Q02846', 'Q03135', 'Q05940', 'Q07075', 'Q12884',
            'Q12933', 'Q13114', 'Q13520', 'Q13740', 'Q13936', 'Q14123',
            'Q14126', 'Q14524', 'Q15628', 'Q16342', 'Q16572', 'Q16625',
            'Q30201', 'Q5VU97', 'Q6NYC1', 'Q7Z7D3', 'Q8NCS7', 'Q8NEC5',
            'Q8WWX8', 'Q92673', 'Q96AP7', 'Q96JQ0', 'Q96SN7', 'Q99712',
            'Q9NP59', 'Q9NQS3', 'Q9NR82', 'Q9NRQ2', 'Q9NSA2', 'Q9NZ08',
            'Q9NZV8', 'Q9P2U7', 'Q9UBH0', 'Q9UKR5', 'Q9ULB1', 'Q9ULH0',
            'Q9UMD9', 'Q9UQD0', 'Q9Y4C0', 'Q9Y5R2', 'Q9Y5Y9', 'Q9Y624',
            'Q9Y6Y9', 'P31997', 'Q5TAT6', 'Q9ULT6', 'O75923',
        },
    ),  # wrong annotations added to exclude
    af.AnnotDef(
        name = 'receptor',
        resource = 'Baccin2019',
        args = {
            'mainclass': 'receptor',
        },
        scope = 'generic',
        # `exclude` defined in `excludes` dict
    ),  # contains also adhesion and some wrong annotations
    af.AnnotDef(
        name = 'receptor',
        resource = 'SignaLink_function',
        args = {
            'function': 'Receptor',
        },
        scope = 'generic',
        exclude = {
            'P48552', 'Q04724', 'P31749', 'P48552', 'P61160',
        },
    ),
    af.AnnotDef(
        name = 'receptor',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Receptor',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'receptor',
        resource = 'Integrins',
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'receptor',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Integrin',
        },
    ),
    # receptor subclasses from HGNC
    af.AnnotDef(
        name = 'immunoglobulin_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Activating leukocyte immunoglobulin like receptors',
        },
    ),
    af.AnnotDef(
        name = 'adiponectin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Adiponectin receptors',
        },
    ),
    af.AnnotDef(
        name = 'adrenalin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Adrenoceptors',
        },
    ),
    af.AnnotDef(
        name = 'angiotensin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Angiotensin receptors',
        },
    ),
    af.AnnotDef(
        name = 'vasopressin_oxytocin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Arginine vasopressin and oxytocin receptors',
        },
    ),
    af.AnnotDef(
        name = 'atypical_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Atypical chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'basigin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Basigin family',
        },
    ),
    af.AnnotDef(
        name = 'bombesin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bombesin receptors',
        },
    ),
    af.AnnotDef(
        name = 'bradykinin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bradykinin receptors',
        },
    ),
    af.AnnotDef(
        name = 'butyrophilin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Butyrophilins',
        },
    ),
    af.AnnotDef(
        name = 'cc_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'C-C motif chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'cx3c_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'C-X-3-C motif chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'cxc_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'C-X-C motif chemokine receptors',
        },
    ),
    af.AnnotDef(
        name = 'celsr_cadherin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'CELSR cadherins',
        },
    ),
    af.AnnotDef(
        name = '5ht_gprotein',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': '5-hydroxytryptamine receptors, G protein-coupled',
        },
    ),
    af.AnnotDef(
        name = '5ht_ionotropic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': '5-hydroxytryptamine receptors, ionotropic',
        },
    ),
    af.AnnotDef(
        name = 'activating_leukocyte_ig_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Activating leukocyte immunoglobulin like receptors',
        },
    ),
    af.AnnotDef(
        name = 'adenosine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Adenosine receptors',
        },
    ),
    af.AnnotDef(
        name = 'calcitonin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Calcitonin receptors',
        },
    ),
    af.AnnotDef(
        name = 'calcium',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Calcium sensing receptors',
        },
    ),
    af.AnnotDef(
        name = 'cannabinoid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cannabinoid receptors',
        },
    ),
    af.AnnotDef(
        name = 'chemerin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chemerin receptor',
        },
    ),
    af.AnnotDef(
        name = 'cholecystokinin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cholecystokinin receptors',
        },
    ),
    af.AnnotDef(
        name = 'muscarinic_cholinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cholinergic receptors muscarinic',
        },
    ),
    af.AnnotDef(
        name = 'nicotinic_cholinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cholinergic receptors nicotinic subunits',
        },
    ),
    af.AnnotDef(
        name = 'collectin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collectins',
        },
    ), # innate immunity receptors for sugar and lipid patterns
    af.AnnotDef(
        name = 'complement_gpcr',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Complement component GPCRs',
        },
    ), # receptors for chemotactic immune signals
    af.AnnotDef(
        name = 'crh',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Corticotropin releasing hormone receptors',
        },
    ),
    af.AnnotDef(
        name = 'dopamine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Dopamine receptors',
        },
    ),
    af.AnnotDef(
        name = 'ephrin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'EPH receptors',
        },
    ),
    af.AnnotDef(
        name = 'endothelin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Endothelin receptors',
        },
    ),
    af.AnnotDef(
        name = 'erbb_rtk',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Erb-b2 receptor tyrosine kinases',
        },
    ),
    af.AnnotDef(
        name = 'f2r',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'F2R receptors',
        },
    ), # GPCRs for thrombin and trypsin
    af.AnnotDef(
        name = 'formyl_peptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Formyl peptide receptors',
        },
    ), # formyl-methionyl peptides are neutrophil chemoattractants
    af.AnnotDef(
        name = 'free_fatty_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Free fatty acid receptors',
        },
    ),  # intestinal short chain fatty acid GPCRs, regulating
        # whole-body energy homeostasis
    af.AnnotDef(
        name = 'bile_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'G protein-coupled bile acid receptor',
        },
    ),  # GPCR for bile acid
    af.AnnotDef(
        name = 'estrogen_gpcr',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'G protein-coupled estrogen receptor',
        },
    ),  # although this receptor is intracellular, there is no reason we
        # shouldn't treat it the same way as plasma membrane receptors
    af.AnnotDef(
        name = 'nuclear_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Nuclear hormone receptors',
        },
    ),  # although these receptors are intracellular, there is no reason we
        # shouldn't treat it the same way as plasma membrane receptors
    af.AnnotDef(
        name = 'gpcr_orphan',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'G protein-coupled receptors, Class A orphans',
                'G protein-coupled receptors, Class C orphans',
            },
        },
    ),  # GPCRs mostly without known ligand, all in the cell membrane
    af.AnnotDef(
        name = 'frizzled',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'G protein-coupled receptors, Class F frizzled',
        },
    ),  # GPCRs for Wnt
    af.AnnotDef(
        name = 'galanin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Galanin receptors',
        },
    ),
    af.AnnotDef(
        name = 'gaba',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Gamma-aminobutyric acid type A receptor subunits',
                'Gamma-aminobutyric acid type B receptor subunits',
            }
        },
    ),
    af.AnnotDef(
        name = 'glucagon',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glucagon receptor family',
        },
    ),
    af.AnnotDef(
        name = 'glutamate_ionotropic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Glutamate ionotropic receptor AMPA type subunits',
                'Glutamate ionotropic receptor NMDA type subunits',
                'Glutamate ionotropic receptor delta type subunits',
                'Glutamate ionotropic receptor kainate type subunits',
            }
        },
    ),
    af.AnnotDef(
        name = 'glutamate_metabotropic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glutamate metabotropic receptors',
        },
    ),
    af.AnnotDef(
        name = 'glycine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glycine receptors',
        },
    ),
    af.AnnotDef(
        name = 'glycoprotein_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glycoprotein hormone receptors',
        },
    ),  # receptors for TSH, FSH and LH
    af.AnnotDef(
        name = 'gonadotropin_releasing_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Gonadotropin releasing hormone receptors',
        },
    ),  # receptors for GnRH
    af.AnnotDef(
        name = 'natriuretic_peptide',
        parent = 'receptor',
        resource = {'P16066', 'P20594'},
    ),
    af.AnnotDef(
        name = 'guanilyn',
        parent = 'receptor',
        resource = {'P25092'},
    ),
    af.AnnotDef(
        name = 'histamine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Histamine receptors',
        },
    ),
    af.AnnotDef(
        name = 'hydroxycarboxylic_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hydroxy-carboxylic acid receptors',
        },
    ),  # receptors for lactate, niacin, etc
    af.AnnotDef(
        name = 'orexin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hypocretin receptors',
        },
    ),
    af.AnnotDef(
        name = 'inhibitory_leukocyte_ig_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Inhibitory leukocyte immunoglobulin like receptors',
        },
    ),
    af.AnnotDef(
        name = 'interferon',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon receptors',
        },
    ),
    af.AnnotDef(
        name = 'interleukin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukin receptors',
        },
    ),
    af.AnnotDef(
        name = 'killer_cell_ig_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Killer cell immunoglobulin like receptors',
        },
    ),  # receptors for HLAs
    af.AnnotDef(
        name = 'killer_cell_lectin_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Killer cell lectin like receptors',
        },
    ),  # receptors mostly for HLAs
    af.AnnotDef(
        name = 'ly6_plaur',
        parent = 'receptor',
        resource = {
            'Q03405', 'Q8IV16',
        }
    ),
    af.AnnotDef(
        name = 'ly6_plaur',
        parent = 'receptor_regulator',
        resource = {
            'Q5SQ64', 'Q8N2G4', 'Q86Y78', 'Q8N6Q3', 'Q16553', 'P0DP58',
            'O94772', 'P13987', 'P0DP57', 'O43653', 'Q8NI32', 'P0C8F1',
        },
    ),  #
    af.AnnotDef(
        name = 'leukotriene',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Leukotriene receptors',
        },
    ),  # receptors for leukotrienes, eicosanoids, N-formyl-met peptides, etc
    af.AnnotDef(
        name = 'pentraxin',
        parent = 'receptor',
        resource = {'Q15818', 'O95502'},
    ),
    af.AnnotDef(
        name = 'short_pentraxin',
        parent = 'secreted_receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Short pentraxins',
        },
    ),  # bind to microbial antigens, DNA and histones
    af.AnnotDef(
        name = 'low_density_lipoprotein',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Low density lipoprotein receptors',
        },
    ),  # receptors for low density lipoproteins
    af.AnnotDef(
        name = 'lysophosphatidic_acid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Lysophosphatidic acid receptors',
        },
    ),
    af.AnnotDef(
        name = 'melanin_concentrating_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Melanin concentrating hormone receptors',
        },
    ),
    af.AnnotDef(
        name = 'melanocortin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Melanocortin receptors',
        },
    ),
    af.AnnotDef(
        name = 'melatonin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Melatonin receptors',
        },
    ),
    af.AnnotDef(
        name = 'membrane_associated_progesterone',
        parent = 'receptor',
        resource = {'O15173', 'O00264'},
    ),  # PGRMC1 is in SER and microsome membrane, but it does not matter
    af.AnnotDef(
        name = 'neuromedin_u',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuromedin U receptors',
        },
    ),
    af.AnnotDef(
        name = 'neuropeptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuropeptide receptors',
        },
    ),
    af.AnnotDef(
        name = 'neurotensin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurotensin receptors',
        },
    ),
    af.AnnotDef(
        name = 'notch',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Notch receptors',
        },
    ),
    af.AnnotDef(
        name = 'olfactory',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Olfactory receptors, family 1',
                'Olfactory receptors, family 2',
                'Olfactory receptors, family 3',
                'Olfactory receptors, family 4',
                'Olfactory receptors, family 5',
                'Olfactory receptors, family 6',
                'Olfactory receptors, family 7',
                'Olfactory receptors, family 8',
                'Olfactory receptors, family 9',
                'Olfactory receptors, family 10',
                'Olfactory receptors, family 11',
                'Olfactory receptors, family 12',
                'Olfactory receptors, family 13',
                'Olfactory receptors, family 14',
                'Olfactory receptors, family 51',
                'Olfactory receptors, family 52',
                'Olfactory receptors, family 56',
            },
        },
    ),  # odorant receptors
    af.AnnotDef(
        name = 'opioid',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Opioid receptors',
        },
    ),
    af.AnnotDef(
        name = 'opsin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Opsin receptors',
        },
    ),
    af.AnnotDef(
        name = 'oxoglutarate',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Oxoglutarate receptor',
        },
    ),
    af.AnnotDef(
        name = 'p2y_purinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'P2Y receptors',
        },
    ),  # receptors for ADP, ATP, UDP, UTP
    af.AnnotDef(
        name = 'p2x_purinergic',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Purinergic receptors P2X',
        },
    ),  # receptors for ATP
    af.AnnotDef(
        name = 'parathyroid_hormone',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Parathyroid hormone receptors',
        },
    ),
    af.AnnotDef(
        name = 'peptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Peptide receptors',
        },
    ),  # receptors for TRH, motilin, apelin, PrRP, QRFP, ghrelin, etc
    af.AnnotDef(
        name = 'platelet_activating_factor',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Platelet activating factor receptor',
        },
    ),
    af.AnnotDef(
        name = 'plexin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Plexins',
        },
    ),  # receptors for semaphorins
    af.AnnotDef(
        name = 'progestin_and_adipoq',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Progestin and adipoQ receptor family',
        },
        exclude = {'Q6TCH7', 'Q8IY49', 'Q15546'},
    ),  # G protein coupled progesterone and ADIPOQ hormone receptors
    af.AnnotDef(
        name = 'prokineticin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Prokineticin receptors',
        },
    ),
    af.AnnotDef(
        name = 'prostaglandin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Prostaglandin receptors',
        },
    ),
    af.AnnotDef(
        name = 'receptor_tyrosine_phosphatase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Protein tyrosine phosphatases receptor type',
        },
        exclude = {'Q92932'},
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'receptor',
        resource = {'P16070', 'Q6UVK1'},
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'receptor_regulator',
        resource = {'O00468'},
    ),
    af.AnnotDef(
        name = 'pseudoautosomal_region',
        parent = 'receptor',
        resource = {'P15509', 'Q86VZ1', 'Q9HC73', 'P26951', 'Q01113'},
    ),
    af.AnnotDef(
        name = 'relt',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'RELT family',
        },
    ),
    af.AnnotDef(
        name = 'gpcr_activity_modifying',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Receptor (G protein-coupled) activity modifying proteins'
            ),
        },
    ),  # receptors for adrenomedullin
    af.AnnotDef(
        name = 'receptor_transporter',
        parent = 'receptor_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Receptor transporter proteins',
        },
    ),  # regulate GPCRs, especially taste and olfactory receptors
    af.AnnotDef(
        name = 'tyrosine_kinase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Receptor tyrosine kinases',
        },
    ),
    af.AnnotDef(
        name = 'relaxin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Relaxin family peptide receptors',
        },
    ),
    af.AnnotDef(
        name = 'repulsive_guidance',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Repulsive guidance molecule family',
        },
    ),  # BMP coreceptors
    af.AnnotDef(
        name = 'slitrk',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'SLIT and NTRK like family',
        },
    ),  # receptors regulating synapse development in CNS
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'receptor',
        resource = {'Q9UEW3', 'Q6ZMJ2', 'P21757', 'P06127', 'Q86VB7'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'secreted_receptor',
        resource = {'Q9UGM3', 'A1L4H1', 'Q86VB7'},
    ),
    af.AnnotDef(
        name = 'plexin',
        parent = 'receptor',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Plexin',
        },
    ),
    af.AnnotDef(
        name = 'ige',
        parent = 'receptor',
        resource = {'Q01362'},
    ),
    af.AnnotDef(
        name = 'scavenger',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Scavenger receptors',
        },
        exclude = {'Q8WTU2', 'Q6AZY7', 'A1L4H1', 'Q14108'}
    ),
    af.AnnotDef(
        name = 'sialic_acid_binding_lectin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sialic acid binding Ig like lectins',
        },
    ),  # not all of them are receptors, some of them might be cell-cell
        # adhesion proteins, depending on the intracellular domains
    af.AnnotDef(
        name = 'somatostatin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Somatostatin receptors',
        },
    ),
    af.AnnotDef(
        name = 'sphingosine_phosphate',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sphingosine 1-phosphate receptors',
        },
    ),
    af.AnnotDef(
        name = 'succinate',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Succinate receptor',
        },
    ),
    af.AnnotDef(
        name = 't_cell',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'T cell receptor alpha locus at 14q11.2',
                'T cell receptor beta locus at 7q34',
                'T cell receptor delta locus at 14q11.2',
                'T cell receptor gamma locus at 7p14',
            },
        },
    ),
    af.AnnotDef(
        name = 'tir_domain',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'TIR domain containing',
        },
        exclude = {'Q8IUC6', 'Q99836', 'Q6SZW1', 'Q86XR7'}
    ),
    af.AnnotDef(
        name = 'tachykinin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tachykinin receptors',
        },
    ),
    af.AnnotDef(
        name = 'taste',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Taste 1 receptors',
                'Taste 2 receptors',
            },
        },
    ),
    af.AnnotDef(
        name = 'toll_like',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Toll like receptors',
        },
    ),
    af.AnnotDef(
        name = 'trace_amin',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Trace amine receptors',
        },
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tumor necrosis factor receptor superfamily',
        },
        exclude = {'O95407', 'O00300'},
    ),
    af.AnnotDef(
        name = 'type1_serine_threonine_kinase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type 1 receptor serine/threonine kinases',
        },
    ),
    af.AnnotDef(
        name = 'type2_serine_threonine_kinase',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type 2 receptor serine/threonine kinases',
        },
    ),
    af.AnnotDef(
        name = 'vasoactive_intestinal_peptide',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Vasoactive intestinal peptide receptor family',
        },
    ),
    af.AnnotDef(
        name = 'vomeronasal',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Vomeronasal receptors',
        },
    ),
    af.AnnotDef(
        name = 'xc_motif_chemokine',
        parent = 'receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'X-C motif chemokine receptors',
        },
    ),
    # subclasses from Almen 2009
    af.AnnotDef(
        name = 'transforming_growth_factor',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Act.TGFB',
        },
    ),
    af.AnnotDef(
        name = 'axl',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Axl',
        },
    ),  # their ligands are in the ECM
    af.AnnotDef(
        name = 'butyrophilin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Butyrophylin',
        },
    ),  # not clear if these are all receptors or some of them are ligands,
        # transporters or regulators of other membrane proteins
    af.AnnotDef(
        name = 'scavenger',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'SCAR',
        },
    ),
    af.AnnotDef(
        name = 'cytokine',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'CytokineR',
        },
    ),
    af.AnnotDef(
        name = 'egf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'EGFR',
        },
    ),
    af.AnnotDef(
        name = 'ephrin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Eph',
        },
    ),
    af.AnnotDef(
        name = 'fgf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'FGFR',
        },
    ),
    af.AnnotDef(
        name = 'fc',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'FcR',
        },
    ),
    af.AnnotDef(
        name = 'frizzled',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Frizzled',
        },
    ),
    af.AnnotDef(
        name = 'gpcr',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'GPCR',
        },
    ),
    af.AnnotDef(
        name = 'glutamate',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Glutamate',
        },
    ),
    af.AnnotDef(
        name = 'ig_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptors',
            'classes': 'IG',
        },
    ),
    af.AnnotDef(
        name = 'il17',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptors',
            'classes': 'IL17',
        },
    ),
    af.AnnotDef(
        name = 'igf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'InsR',
        },
    ),
    af.AnnotDef(
        name = 'killer_cell_ig_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'KIR',
        },
    ),
    af.AnnotDef(
        name = 'kinase',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'mainclass': 'Receptors',
            'classes': {
                'Kinase',
                'KInase' # this is a typo
            },
        },
    ),  # all receptors with kinase activity,
        # I think it's useful to have such category
    af.AnnotDef(
        name = 'low_density_lipoprotein',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'LDLR',
        },
    ),
    af.AnnotDef(
        name = 'leukocyte_ig_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'LILR',
        },
    ),
    af.AnnotDef(
        name = 'mannose',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'MacrophageMannoseR',
        },
    ),
    af.AnnotDef(
        name = 'netrin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'NetrinR',
        },
    ),
    af.AnnotDef(
        name = 'neuropilin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Neuropilin',
        },
    ),  # receptor for semaphorins, VEGF and PLGF
    af.AnnotDef(
        name = 'notch',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Notch',
        },
    ),
    af.AnnotDef(
        name = 'olfactory',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': {'Olf', 'PutativeOlfR'},
        },
    ),
    af.AnnotDef(
        name = 'cd1',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherCD1',
        },
    ),
    af.AnnotDef(
        name = 'cd300',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherCD300',
        },
    ),
    af.AnnotDef(
        name = 'natural_cytotoxicity_triggering',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherNCR',
        },
    ),
    af.AnnotDef(
        name = 'poliovirus',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherPVR',
        },
    ),
    af.AnnotDef(
        name = 'roundabout',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherROBO',
        },
    ),
    af.AnnotDef(
        name = 'triggering',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'OtherTREM',
        },
        exclude = {'Q6UXN2'}
    ),
    af.AnnotDef(
        name = 'pdgf',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'PDGFR',
        },
    ),
    af.AnnotDef(
        name = 'patched',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Patched',
        },
    ),
    af.AnnotDef(
        name = 'plexin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Plexin',
        },
    ),
    af.AnnotDef(
        name = 'receptor_tyrosine_phosphatase',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'ReceptorTypePhosphatases',
        },
    ),  # not sure all these are receptors
    af.AnnotDef(
        name = 'rhodopsin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Rhodopsin',
        },
    ),  # GPCRs for various ligands, hormones, neurotransmitters, etc
    af.AnnotDef(
        name = 'secretin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Secretin',
        },
    ),  # GPCRs for various hormones
    af.AnnotDef(
        name = 'syndecan',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Syndecan',
        },
    ),  # heparan sulphate carrying cell surface proteins
        # transferring signals to the cytoskeleton, regulating cell shape
    af.AnnotDef(
        name = 'taste',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TAS2R',
        },
    ),
    af.AnnotDef(
        name = 't_cell',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TCR',
        },
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TNFNGF',
        },
    ),  # maybe this is a broader superfamily
    af.AnnotDef(
        name = 'toll_like',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'TOLL',
        },
    ),  # receptors for various pathogen patterns e.g. LPS
    af.AnnotDef(
        name = 'teneurin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Teneurin',
        },
    ),  # functioning in neural development and maybe elsewhere;
        # these are definitely receptors, but also ligands and
        # maybe cell-cell adhesion molecules
    af.AnnotDef(
        name = 'transferrin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Transferrin',
        },
    ),
    af.AnnotDef(
        name = 'type1_ig_like_cytokine',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Type1',
        },
    ),  # mostly interleukin receptors
    af.AnnotDef(
        name = 'type2_ig_like_cytokine',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'Type2',
        },
    ),  # mostly interferon and interleukin receptors
    af.AnnotDef(
        name = 'vomeronasal',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'V1R',
        },
    ),  # pheromone receptors
    af.AnnotDef(
        name = 'gaba_ach',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'cys-loop',
        },
    ),  # receptors for ACh and GABA
    af.AnnotDef(
        name = 'neurotrophin',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'neutrophin',
        },
    ),
    # specific classes from UniProt keywords
    af.AnnotDef(
        name = 'gpcr',
        parent = 'receptor',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'G-protein coupled receptor',
        },
    ),
    af.AnnotDef(
        name = 'retinal_guanylyl_cyclase',
        parent = 'receptor',
        resource = 'Almen2009',
        args = {
            'classes': 'RGC',
        },
    ),  # receptors for various compounds, e.g. natriuretic peptide
        # and E.coli enterotoxin
    af.AnnotDef(
        name = 'receptor_transporter',
        parent = 'receptor_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'RTP',
        },
    ),  # regulate other receptors,
        # especially their trafficking to the plasma membrane
    af.AnnotDef(
        name = 'lectin',
        parent = 'receptor',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Lectin',
        },
        limit = 'plasma_membrane_transmembrane',
        avoid = af.AnnotDef(
            name = 'ecm_located',
            resource = 'UniProt_location',
            args = {
                'location': 'Extracellular matrix',
            },
        ),
    ),


    # secreted receptors
    af.AnnotDef(
        name = 'secreted_receptor',
        source = 'composite',
        scope = 'generic',
        resource = '~secreted_receptor',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'ly6_plur',
        parent = 'secreted_receptor',
        resource = {'Q6UX82', 'P55000', 'P13987'},
    ),
    af.AnnotDef(
        name = 'pentraxin',
        parent = 'secreted_receptor',
        resource = {'Q96A99', 'P26022', 'P47972'},
    ),
    af.AnnotDef(
        name = 'peptidoglycan',
        parent = 'secreted_receptor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Peptidoglycan recognition proteins',
        },
    ),  # apart from peptide recognition they have anti-microbial activity
        # either enzymatically or other ways
    af.AnnotDef(
        name = 'scavenger',
        parent = 'secreted_receptor',
        resource = {'Q8WTU2', 'A1L4H1', 'Q86VB7'},
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'secreted_receptor',
        resource = {'O9540', 'O0030'},
    ),
    af.AnnotDef(
        name = 'lectin',
        parent = 'secreted_receptor',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Lectin',
        },
        limit = 'secreted',
        avoid = af.AnnotDef(
            name = 'ecm',
            resource = 'UniProt_location',
            args = {
                'location': 'Extracellular matrix',
            },
        ),
    ),
    af.AnnotDef(
        name = 'galectin',
        parent = 'secreted_receptor',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Galectin',
        },
    ),
    af.AnnotDef(
        name = 'secreted_receptor',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Receptor',
        },
        limit = 'secreted',
    ),

    # receptor regulators
    af.AnnotDef(
        name = 'receptor_regulator',
        source = 'composite',
        scope = 'generic',
        resource = '~receptor_regulator',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'activating_cofactor',
        parent = 'receptor_regulator',
        resource = 'CellChatDB',
        args = {
            'role': 'co_A_receptor',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'inhibitory_cofactor',
        parent = 'receptor_regulator',
        resource = 'CellChatDB',
        args = {
            'role': 'co_I_receptor',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ms4',
        parent = 'receptor_regulator',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'ms4',
                    resource = 'HGNC',
                    args = {
                        'mainclass': 'Membrane spanning 4-domains',
                    },
                ),
                'plasma_membrane',
            ),
            op = set.intersection,
        ),
    ),
    af.AnnotDef(
        name = 'receptor_activity_modifying',
        parent = 'receptor_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'RAMP',
        },
    ),  # affect other receptors such as CALCRL
    af.AnnotDef(
        name = 'glypican',
        parent = 'receptor_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glypicans',
        },
    ),
    af.AnnotDef(
        name = 'glypican',
        parent = 'receptor_regulator',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Glypican',
        },
    ),
    af.AnnotDef(
        name = 'syndecan',
        parent = 'receptor_regulator',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Syndecan',
        },
    ),

    # ECM
    af.AnnotDef(
        name = 'ecm',
        resource = '~ecm',
        scope = 'generic',
        source = 'composite',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'ecm',
        resource = 'CellChatDB',
        args = {
            'role': 'ligand',
            'category': 'ECM-Receptor',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ecm',
        resource = 'Cellinker',
        args = {
            'location': 'ECM',
        },
        scope = 'generic',
    ),
    # ecm: UniProt specific categories
    af.AnnotDef(
        name = 'lectin',
        parent = 'ecm',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Lectin',
        },
        limit = af.AnnotDef(
            name = 'ecm',
            resource = 'UniProt_location',
            args = {
                'location': 'Extracellular matrix',
            },
        ),
    ),
    af.AnnotDef(
        name = 'collagen',
        parent = 'ecm',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Lectin',
        },
        limit = af.AnnotDef(
            name = 'ecm',
            resource = 'UniProt_location',
            args = {
                'location': 'Extracellular matrix',
            },
        ),
    ),
    # ecm: Matrisome specific categories
    af.AnnotDef(
        name = 'collagen',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subclass': 'Collagens',
        },
    ),
    af.AnnotDef(
        name = 'glycoprotein',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subclass': 'ECM Glycoproteins',
        },
        exclude = {
            'Q8TC99', 'Q502W6', 'Q96HD1', 'P55081', 'Q8IUX7', 'Q92832',
            'Q6UXH1', 'Q5JSJ4', 'Q96SY0', 'Q08431', 'Q9Y215', 'O95389',
            'A1KZ92', 'Q8WWZ8',
        },
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subclass': 'Proteoglycans',
        },
        exclude = {'P10124', 'Q9Y2Y8'},
    ),
    af.AnnotDef(
        name = 'mucin',
        parent = 'ecm',
        resource = {
            'Q5SZK8', 'Q6UVK1', 'P35247', 'Q99102', 'Q9Y625', 'P09382',
            'P0C091', 'Q9ULC0', 'Q8N387', 'Q5H8C1', 'E2RYF6', 'Q02817',
            'Q8TAX7', 'P98088', 'Q9HC84', 'Q6W4X9', 'Q7Z5P9', 'Q9UKN1',
            'Q685J3', 'Q9H3R2', 'Q8WXI7', 'P15941', 'Q8N387', 'Q8N307',
            'Q02505', 'Q5SSG8',
        },
    ),  # these are mostly mucins, selected by Leila from the
        # Matrisome "ECM-affiliated" category; however this category
        # contains 6x more proteins and later would be good to review
        # these again and include in further categories
    af.AnnotDef(
        name = 'basement_membrane',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Basement Membrane',
        },
    ),
    af.AnnotDef(
        name = 'fibril_associated_collagen_with_interrupted_triple_helices',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'FACIT',
        },
    ),
    af.AnnotDef(
        name = 'fibulin',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Fibulin',
        },
    ),
    af.AnnotDef(
        name = 'hemostasis_related',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Hemostasis',
        },
    ),
    af.AnnotDef(
        name = 'laminin',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Laminin, Basement Membrane',
        },
    ),
    af.AnnotDef(
        name = 'mucin',
        parent = 'ecm',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Mucin',
        },
    ),
    af.AnnotDef(
        name = 'ecm',
        scope = 'generic',
        resource = 'MatrixDB',
        args = {
            'mainclass': 'ecm',
        },
        exclude = {
            'Q12959', 'Q96JX3', 'Q96SY0', 'Q9UQE7',
            'O95389', 'Q9NZU0', 'P48745', 'O75445',
            'O75900', 'P05230', 'Q8TC99', 'O43155',
            'Q03001', 'Q5JSJ4', 'Q53GQ0', 'Q9Y5L3',
            'Q14129', 'A1KZ92', 'Q9Y2Y8', 'Q9UBB9',
            'Q9NZU1', 'P27797', 'Q92896', 'Q96HD1',
            'Q8IVL5', 'P22303', 'P23229', 'Q96RT1',
            'Q7RTW8', 'Q502W6', 'P55081', 'Q8IUX7',
            'Q92832', 'Q6UXH1', 'Q08431', 'Q9Y215',
            'Q03167', 'P03950', 'P10124', 'P52803',
            'Q96SY0', 'Q5KU26', 'Q86UN3', 'Q9Y2Y8',
            'O15455', 'Q96HD1', 'Q86UN2', 'Q6EMK4',
            'Q502W6', 'Q6AZY7', 'Q92832', 'Q6UXH1',
            'Q7L0X0', 'P10124',
        },
        avoid = ('ligand', 'secreted_enzyme', 'receptor'),
    ),  # some potentially wrong elements such as ligands, with removal
        # of the groups above it looks more or less fine
    af.AnnotDef(
        name = 'ecm',
        scope = 'generic',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'ecm structure',
        },
    ),
    af.AnnotDef(
        name = 'ecm',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'location': 'extracellular matrix',
        },
        exclude = {
            'P55081', 'P03950', 'Q9UBV4', 'O14904', 'O00744', 'Q9H1J7',
            'O00755', 'Q9NPA2', 'P56706', 'P09544', 'Q93097', 'P10745',
            'Q9GZT5', 'O14905', 'O96014', 'Q9NRE1', 'P56703', 'P01009',
            'P04628', 'P09238', 'Q9Y6F9', 'P08254', 'P56704', 'P39900',
            'P01137', 'Q5K4E3', 'P03956', 'Q9H239', 'P22894', 'Q92626',
            'P56705', 'P41221', 'Q93098', 'Q9H1J5', 'Q6UY14', 'Q66K79',
            'P55081', 'P08253', 'P09237', 'P12644', 'P13497', 'P45452',
            'Q8N2S1', 'P35625', 'P14780', 'Q8N6G6', 'Q8IUX8', 'Q9ULZ9',
            'P51512',
        },
    ),
    af.AnnotDef(
        name = 'basement membrane',
        parent = 'ecm',
        resource = 'Ramilowski_location',
        args = {
            'location': 'basement membrane',
        },
    ),

    af.AnnotDef(
        name = 'ecm',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Extracellular matrix',
        },
    ),
    af.AnnotDef(
        name = 'ecm',
        resource = 'CellCellInteractions',
        scope = 'generic',
        args = {
            'mainclass': 'ECM',
        },
        exclude = {
            'Q12959', 'Q01196', 'Q6PIL6', 'Q9NS61',
            'Q96JX3', 'Q9UQE7', 'O95389', 'Q9NZU0',
            'P48745', 'O75445', 'O75197', 'O43294',
            'O75900', 'Q8N7U6', 'P05230', 'P16144',
            'O43155', 'Q05586', 'P00441', 'Q03001',
            'P07900', 'P26012', 'Q53GQ0', 'Q9Y5L3',
            'P14618', 'P37840', 'Q4KMG0', 'Q14129',
            'P16150', 'Q14651', 'Q9UBB9', 'Q9NZU1',
            'P07339', 'P27797', 'Q92896', 'P08311',
            'Q9NZU5', 'P16591', 'Q8IVL5', 'P22303',
            'Q9ULV1', 'P23229', 'P05556', 'Q8IVL1',
            'P21802', 'Q96RT1', 'Q7RTW8', 'Q9NZI2',
            'P55081', 'P29122', 'P52272', 'Q9Y215',
            'P14625', 'Q03167', 'Q9Y2W7', 'P06756',
            'Q01638', 'P03950', 'Q9H1J7',
        },
        avoid = ('ligand', 'receptor', 'secreted_enzyme'),
    ),  # more or less correct, but includes enzymes and matrix adhesion
        # these we excluded
        # TODO: the rest can be added to ligand and secreted enzyme
        # categories
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'collagen_proteoglycan',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collagen proteoglycans',
        },
    ),
    af.AnnotDef(
        name = 'collagen',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collagens',
        },
    ),
    af.AnnotDef(
        name = 'emi',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'EMI domain containing',
        },
    ),  # this could be also cell-matrix adhesion, although these proteins
        # are not in the cell membrane but all secreted
    af.AnnotDef(
        name = 'fibrillin',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibrillins',
        },
    ),
    af.AnnotDef(
        name = 'laminin',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Laminin subunits',
        },
    ),
    af.AnnotDef(
        name = 'fibulin',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibulins',
        },
    ),  # parts of ECM, especially elastic fibers, one of them is a ligand
        # for EGFR (but still an EVM protein at the same time)
    af.AnnotDef(
        name = 'hyalectan_proteoglycan',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hyalectan proteoglycans',
        },
    ),
    af.AnnotDef(
        name = 'matrilin',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Matrilins',
        },
    ),  # cartilage ECM
    af.AnnotDef(
        name = 'mucin',
        parent = 'ecm',
        resource ='HGNC',
        args = {
            'mainclass': 'Mucins',
        },
        limit = 'secreted',
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'ecm',
        resource = {'O00468', 'P98160'},
    ),
    af.AnnotDef(
        name = 'sibling',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'SIBLING family',
        },
    ),
    af.AnnotDef(
        name = 'sparc_ecm_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'SPARC family',
        },
    ),  # act either on ligands or ECM or both
    af.AnnotDef(
        name = 'small_leucine_rich_repeat_proteoglycan',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Small leucine rich repeat proteoglycans',
        },
    ),
    af.AnnotDef(
        name = 'zona_pellucida_glycoprotein',
        parent = 'ecm',
        resource = 'HGNC',
        args = {
            'mainclass': 'Zona pellucida glycoproteins',
        },
    ),  # ECM of the zona pellucida (zone surrounding the oocyte)

    # ECM regulators
    af.AnnotDef(
        name = 'ecm_regulator',
        resource = af.AnnotOp(
            annots = '~ecm_regulator',
            op = set.union,
        ),
        receiver = False,
        transmitter = True,
        scope = 'generic',
        source = 'composite',
    ),
    af.AnnotDef(
        name = 'ecm_regulator',
        scope = 'generic',
        resource = 'Matrisome',
        args = {
            'subclass': 'ECM Regulators',
        },
        exclude = {
            'O00469', 'O15460', 'O43548', 'O60911', 'O75063', 'O75635',
            'O95932', 'P00488', 'P01040', 'P04080', 'P07339', 'P09668',
            'P10619', 'P13674', 'P14091', 'P20848', 'P29508', 'P35237',
            'P43234', 'P48594', 'P48595', 'P50452', 'P50453', 'P50454',
            'P53634', 'P56202', 'Q08188', 'Q6HA08', 'Q6YHK3', 'Q7Z4N8',
            'Q86WD7', 'Q8IVL5', 'Q8IVL6', 'Q8NBH2', 'Q96IV0', 'Q96KS0',
            'Q96P15', 'Q9GZT9', 'Q9H6Z9', 'Q9NXG6', 'Q9UBR2', 'Q9UBX1',
            'Q9UIV8', 'Q9UKF2',
        },
    ),  # mostly secreted enzymes acting on ECM components,
        # proteases and protease inhibitors

    # ligand
    af.AnnotDef(
        name = 'ligand',
        resource = '~ligand',
        receiver = False,
        transmitter = True,
        scope = 'generic',
        source = 'composite',
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'talklr',
        args = {
            'role': 'ligand',
            'putative': False,
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'Cellinker',
        args = {
            'role': 'ligand',
            'type': {
                'Cytokine-cytokine receptor interaction',
                'Secreted protein to receptor interaction',
            },
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'scConnect',
        args = {
            'role': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'connectomeDB2020',
        args = {
            'role': 'ligand',
            'location': 'secreted',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'Matrisome',
        args = {
            'subclass': 'Secreted Factors',
        },
        exclude = {'Q14512', 'P51610'},
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'CellCall',
        args = {
            'role': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'cytokine',
        parent = 'ligand',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Cytokine',
        },
        limit = 'secreted',
    ),
    af.AnnotDef(
        name = 'growth_factor',
        parent = 'ligand',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Growth factor',
        },
        exclude = {'Q6ZN28', 'P26441', 'Q9Y3E1'},
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'iTALK',
        scope = 'generic',
        args = {
            'mainclass': 'ligand',
        },
        exclude = {
            'P39019', 'P07237', 'P19823', 'P17752', 'Q86UR5', 'P09211',
            'P02741', 'P0DP24', 'P05543', 'O60494', 'P05067', 'P02679',
            'P05155', 'Q96DA0', 'P84077', 'P07900', 'P78536', 'P80188',
            'O95467', 'O75596', 'P35555', 'P0DP25', 'Q13444', 'Q9NRV9',
            'P35354', 'P07288', 'O75077', 'Q2MV58', 'P08670', 'P00734',
            'P00748', 'Q5JWF2', 'O95711', 'P10153', 'Q8N474', 'P23515',
            'P00797', 'P07602', 'Q9Y2I2', 'P16520', 'P43489', 'P00742',
            'P09417', 'Q4VX76', 'P11226', 'Q8WWY8', 'P43490', 'Q9NUP9',
            'O00592', 'P07093', 'O15496', 'P00995', 'P52945', 'Q8TAX7',
            'P12830', 'P02765', 'Q10588', 'O95428', 'Q16613', 'P05121',
            'Q9NSG2', 'P04278', 'Q9UIW2', 'P04196', 'Q8IZJ3', 'Q9NQC3',
            'P84996', 'Q9Y566', 'P19113', 'P07225', 'O14672', 'P08603',
            'O60469', 'P01112', 'O94887', 'P42127', 'Q9BQ66', 'P09429',
            'P09429', 'P03973', 'P06454', 'O95084', 'P19835', 'Q9XRX5',
            'P02671', 'P00749', 'Q5SR53', 'P49768', 'P08861', 'Q13442',
            'P08571', 'P42081', 'P51654', 'P22692', 'Q13352', 'Q9Y215',
            'P04070', 'P01009', 'P01031', 'Q6UWX4', 'P61769', 'P20062',
            'P49913', 'Q9UQ26', 'P01023', 'P00750', 'P02746', 'O43184',
            'Q15165', 'P02675', 'P00488', 'P02753', 'P01033', 'P08174',
            'Q99259', 'P13385', 'P10144', 'P14618', 'Q8N2X6', 'P61626',
            'P54577', 'P00740', 'P00451', 'Q9BX66', 'Q6NW40', 'P0C0L4',
            'Q02817', 'Q6V0I7', 'P01024', 'P01008', 'P21815', 'P03951',
            'O00587', 'Q07326', 'P31025', 'P14625', 'P02768', 'P48651',
            'P12821', 'O14594', 'P10646', 'P63092', 'Q92819', 'P02788',
            'Q9UKQ2', 'P17752', 'P19823', 'P07237', 'P39019', 'P50552',
            'Q86UR5', 'Q92956', 'P55058', 'O95274', 'P80303', 'Q5T5A4',
            'P62987', 'Q9NZR2', 'P21980', 'P08709', 'P16035', 'O43278',
            'P20292', 'P08582', 'Q13477', 'P25063', 'Q7Z6A9', 'P0DP23',
            'Q9NV23',
        }
    ),  # locations are 90% correct (secreted or cell surface), but includes
        # some ECM, enzyme, regulator, etc proteins
    af.AnnotDef(
        name = 'growth_factor',
        parent = 'ligand',
        resource = 'iTALK',
        args = {
            'mainclass': 'ligand',
            'subclass': 'growth factor',
        },
    ),  # looks good
    af.AnnotDef(
        name = 'cytokine',
        parent = 'ligand',
        resource = 'iTALK',
        args = {
            'mainclass': 'ligand',
            'subclass': 'cytokine',
        },
    ),  # looks good
    af.AnnotDef(
        name = 'ccn_family',
        parent = 'ligand',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'CCN Family',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'CellCellInteractions',
        args = {
            'mainclass': 'Ligand',
        },
        scope = 'generic',
        exclude = {
            'Q6UWW8', 'P15531', 'P08236', 'Q9UHG3', 'Q9H9H4', 'O43852',
            'Q92896', 'P55789', 'P22392', 'Q96A49', 'Q8WZ79', 'Q9BS26',
            'O95236', 'Q9UJU6', 'Q8NHP8', 'P35475', 'P12109',
        },
        enabled = False,
    ),  # includes both secreted and cell surface,
        # also enzymes and regulators
    af.AnnotDef(
        name = 'ligand',
        resource = 'EMBRACE',
        args = {
            'mainclass': 'ligand',
        },
        scope = 'generic',
        exclude = {
            'P35354', 'P14618', 'Q4VX76', 'Q2MV58', 'P84077', 'P08709',
            'P12110', 'P06858', 'P11150', 'P04899', 'P02745', 'P20849',
            'P39060', 'P27658', 'P45452', 'P98160', 'Q05707', 'P12107',
            'P02452', 'P08253', 'Q02388', 'P39900', 'P08123', 'Q96CG8',
            'P21941', 'P35625', 'P12109', 'P20908', 'P02461', 'O00587',
            'O14594', 'O14672', 'O15496', 'O43184', 'O43278', 'O60469',
            'O60494', 'O75077', 'O94887', 'O95084', 'O95467', 'O95711',
            'P00451', 'P00488', 'P00734', 'P00740', 'P00742', 'P00749',
            'P00750', 'P00797', 'P00995', 'P01008', 'P01023', 'P01024',
            'P01031', 'P01033', 'P01112', 'P02746', 'P02768', 'P02788',
            'P03973', 'P04278', 'P05067', 'P05155', 'Q9Y566', 'P07093',
            'P07225', 'P07237', 'P07602', 'P07900', 'P08571', 'P08670',
            'P08709', 'P09417', 'P10646', 'P12821', 'P12830', 'P14618',
            'P16520', 'P19823', 'P20062', 'P21815', 'P21980', 'P22692',
            'P23515', 'P35354', 'P35555', 'P39019', 'P42127', 'P43490',
            'P48651', 'P49768', 'P49913', 'P51654', 'P55058', 'P61769',
            'P63092', 'P78536', 'P80188', 'P84077', 'P84996', 'Q07326',
            'Q10588', 'Q13352', 'Q13442', 'Q13444', 'Q15165', 'Q2MV58',
            'Q4VX76', 'Q5JWF2', 'Q6NW40', 'Q6UWX4', 'Q6V0I7', 'Q7Z6A9',
            'Q86UR5', 'Q8N474', 'Q8WWY8', 'Q92819', 'Q9BX66', 'Q9NQC3',
            'Q9NUP9', 'Q9NV23', 'Q9UQ26', 'Q9Y2I2',
        },
    ),  # inlcudes secreted enzymes, some ECM proteins which we exclude
    af.AnnotDef(
        name = 'ligand',
        resource = '~ligand~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'CellPhoneDB',
        scope = 'generic',
        args = {
            'secreted': bool,
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'GO_Intercell',
        scope = 'generic',
        args = {
            'mainclass': 'ligands',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'HPMR',
        scope = 'generic',
        args = {
            'role': 'Ligand',
        },
        # `exclude` defined in the `excludes` dict
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'ICELLNET',
        args = {
            'role': 'ligand',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'CellChatDB',
        args = {
            'role': 'ligand',
            'category': 'Secreted Signaling',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'CellTalkDB',
        args = {
            'role': 'ligand',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'Ramilowski2015',
        args = {
            'mainclass': 'ligand',
        },
        exclude = {
            'O00587', 'O00592', 'O14672', 'O15496', 'O43184', 'O43278',
            'O60494', 'O75077', 'O95084', 'O95467', 'O95711', 'P00451',
            'P00488', 'P00734', 'P00740', 'P00742', 'P00748', 'P00749',
            'P00750', 'P00797', 'P00995', 'P01008', 'P01009', 'P01023',
            'P01024', 'P01031', 'P01033', 'P01112', 'P02452', 'P02461',
            'P02671', 'P02675', 'P02679', 'P02741', 'P02745', 'P02746',
            'P02753', 'P02765', 'P02768', 'P02788', 'P03951', 'P04070',
            'P04196', 'P04278', 'P04899', 'P05067', 'P05121', 'P05155',
            'P05543', 'P06454', 'P04003', 'P06858', 'P07093', 'P07225',
            'P07288', 'P07602', 'P07900', 'P08123', 'P08174', 'P08253',
            'P08571', 'P08603', 'P08670', 'P08709', 'P09211', 'P09429',
            'P0C0L4', 'P10144', 'P10153', 'P10646', 'P11150', 'P11226',
            'P12107', 'P12821', 'P12830', 'P13385', 'P16035', 'P16520',
            'P19823', 'P19835', 'P20292', 'P20849', 'P20908', 'P21941',
            'P21980', 'P22692', 'P23515', 'P25063', 'P27658', 'P31025',
            'P35354', 'P35555', 'P35625', 'P39019', 'P39060', 'P39900',
            'P42127', 'P43490', 'P45452', 'P48651', 'P49768', 'P49913',
            'P50552', 'P51654', 'P54577', 'P55058', 'P61769', 'P63092',
            'P78536', 'P80188', 'P80303', 'P84077', 'P84996', 'P98160',
            'Q02817', 'Q05707', 'Q07326', 'Q10588', 'Q13352', 'Q13442',
            'Q13444', 'Q13477', 'Q15165', 'Q16613', 'Q2MV58', 'Q4VX76',
            'Q5JWF2', 'Q6NW40', 'Q6V0I7', 'Q7Z6A9', 'Q86UR5', 'Q8IZJ3',
            'Q8N2X6', 'Q8N474', 'Q8TAX7', 'Q8WWY8', 'Q92819', 'Q96CG8',
            'Q96DA0', 'Q9BQ66', 'Q9BX66', 'Q9NQC3', 'Q9NRV9', 'Q9NUP9',
            'Q9NZR2', 'Q9UIW2', 'Q9UKQ2', 'Q9UQ26', 'Q9Y215', 'Q9Y2I2',
            'Q9Y566', 'P05997', 'P08572', 'P02458', 'Q8IWL1', 'P29400',
            'P25940', 'Q13443', 'Q14055', 'P03956', 'Q14050', 'Q99965',
            'Q14766', 'P53420',
        },
    ),  # many non-ligand proteins excluded
    af.AnnotDef(
        name = 'ligand',
        scope = 'generic',
        resource = 'Kirouac2010',
        args = {
            'mainclass': 'ligand',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'Guide2Pharma',
        scope = 'generic',
        args = {
            'mainclass': 'ligand',
        },
        exclude = {
            'P54750', 'Q14123', 'P09917', 'Q02846', 'Q01064',
        },
    ),  # great collection of ligands
    af.AnnotDef(
        name = 'ligand',
        resource = af.AnnotOp(
            annots = '~ligand~DGIdb',
            op = set.union,
        ),
        scope = 'generic',
        resource_name = 'DGIdb',
    ),
    af.AnnotDef(
        name = 'growth_factor',
        parent = 'ligand',
        resource = 'DGIdb',
        args = {
            'category': 'GROWTH FACTOR',
        },
        exclude = {
            'P0C7T3', 'P55789', 'P13385', 'P01033', 'Q9UIW2', 'P00734',
            'O14649', 'O75534', 'Q8NGH5', 'Q8NGH8', 'P21917', 'P14416',
            'P26441', 'Q8NGS4', 'Q9H5N1', 'Q8NG75', 'P35462',
        },
    ),  # good, apart from a few surprising exceptions
    af.AnnotDef(
        name = 'hormone',
        parent = 'ligand',
        resource = 'DGIdb',
        args = {
            'category': 'HORMONE ACTIVITY',
        },
        exclude = {
            'P21917', 'P14416', 'P0C7T3', 'O75534', 'Q8NGH8', 'P35555',
            'Q8NGH5', 'P35462', 'Q8NGS4', 'O95751', 'P26441', 'P13385',
            'Q9H5N1', 'Q8NG75', 'Q9UIW2',
        },
    ),  # good, apart from a few surprising exceptions
    af.AnnotDef(
        name = 'ligand',
        resource = 'LRdb',
        scope = 'generic',
        args = {
            'role': 'ligand',
            'references': bool,
        },
        exclude = {
            'O00468', 'O00587', 'O00592', 'O14672', 'O15496', 'O43184',
            'O43278', 'O43914', 'O60494', 'O75077', 'O95084', 'O95467',
            'O95711', 'O95754', 'P00451', 'P00488', 'P00734', 'P00740',
            'P00742', 'P00748', 'P00749', 'P00750', 'P00797', 'P00813',
            'P00995', 'P01008', 'P01009', 'P01023', 'P01024', 'P01031',
            'P01033', 'P01112', 'P02452', 'P02458', 'P02461', 'P02462',
            'P02654', 'P02671', 'P02675', 'P02679', 'P02741', 'P02745',
            'P02746', 'P02753', 'P02765', 'P02768', 'P02788', 'P03951',
            'P03956', 'P04003', 'P04070', 'P04196', 'P04278', 'P04899',
            'P05067', 'P05107', 'P05121', 'P05155', 'P05543', 'P05997',
            'P06454', 'P06734', 'P06858', 'P07093', 'P07225', 'P07288',
            'P07602', 'P07900', 'P08123', 'P08174', 'P08253', 'P08571',
            'P08572', 'P08603', 'P08670', 'P08709', 'P09211', 'P09237',
            'P09429', 'P0C0L4', 'P0CG37', 'P0DMV8', 'P0DP23', 'P10144',
            'P10153', 'P10646', 'P11150', 'P11226', 'P12107', 'P12821',
            'P12830', 'P13385', 'P13591', 'P13688', 'P15151', 'P15907',
            'P16035', 'P16520', 'P16671', 'P19823', 'P19835', 'P20273',
            'P20849', 'P20908', 'P21810', 'P21941', 'P21980', 'P22692',
            'P22897', 'P23515', 'P25063', 'P25940', 'P26441', 'P27658',
            'P29400', 'P30533', 'P31025', 'P35354', 'P35555', 'P35613',
            'P35625', 'P39019', 'P39060', 'P39900', 'P42081', 'P42127',
            'P43405', 'P43490', 'P45452', 'P48651', 'P49768', 'P49913',
            'P50552', 'P51654', 'P53420', 'P55058', 'P61769', 'P63092',
            'P78536', 'P80188', 'P80303', 'P84077', 'P84996', 'P98160',
            'Q01469', 'Q02817', 'Q05707', 'Q07326', 'Q08722', 'Q10588',
            'Q12918', 'Q13158', 'Q13352', 'Q13361', 'Q13442', 'Q13443',
            'Q13444', 'Q13477', 'Q14031', 'Q14050', 'Q14055', 'Q14242',
            'Q14766', 'Q15165', 'Q15223', 'Q16613', 'Q2MV58', 'Q4VX76',
            'Q5JWF2', 'Q6NW40', 'Q6V0I7', 'Q7Z6A9', 'Q86UR5', 'Q8IWL1',
            'Q8IZJ3', 'Q8N474', 'Q8NHJ6', 'Q8TAX7', 'Q8WWY8', 'Q92819',
            'Q92854', 'Q96CG8', 'Q96DA0', 'Q99965', 'Q9BQ66', 'Q9BX66',
            'Q9C0C4', 'Q9H2E6', 'Q9H3S1', 'Q9NQC3', 'Q9NRV9', 'Q9NTN9',
            'Q9NUP9', 'Q9NWZ3', 'Q9NZR2', 'Q9UIW2', 'Q9UKQ2', 'Q9UQ26',
            'Q9Y215', 'Q9Y2I2', 'Q9Y566', 'Q9Y625',
        },
    ),
    af.AnnotDef(
        name = 'ligand',
        resource = 'Baccin2019',
        scope = 'generic',
        args = {
            'mainclass': 'ligand',
            'location': {'secreted', 'both', 'ecm'},
        },
        # `exclude` defined in `excludes` dict
    ),  # quite some wrong annotations, also many receptors
    af.AnnotDef(
        name = 'ligand',
        resource = 'SignaLink_function',
        scope = 'generic',
        args = {
            'function': 'Ligand',
        },
        exclude = {
            'Q08334', 'Q8NFT8', 'P26441', 'O00468', 'Q8IZL2', 'P18564',
            'Q13361', 'P26012',
        },
    ),
    # specific ligand classes from HGNC
    af.AnnotDef(
        name = 'angiopoietin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Angiopoietin like family',
        },
    ),
    af.AnnotDef(
        name = 'basigin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Basigin family',
        },
    ),
    af.AnnotDef(
        name = 'bmp',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bone morphogenetic proteins',
        },
    ),
    af.AnnotDef(
        name = 'c1q',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'C1q and TNF related',
        },
    ),
    af.AnnotDef(
        name = 'ccn',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cellular communication network factors',
        },
    ),
    af.AnnotDef(
        name = 'endogenous',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Endogenous ligands',
        },
    ), # a very few among these are actually not secreted
    af.AnnotDef(
        name = 'chemokine',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chemokine ligands',
        },
    ),
    af.AnnotDef(
        name = 'chordin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chordin family',
        },
    ), # BMP antagonist ligands
    af.AnnotDef(
        name = 'cysteine_rich_bmp_regulator',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cysteine rich transmembrane BMP regulators',
        },
    ), # BMP agonist and antagonist ligands
    af.AnnotDef(
        name = 'dan',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'DAN family',
        },
    ), # TGF & BMP signaling agonists and antagonists
    af.AnnotDef(
        name = 'fgf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibroblast growth factor family',
        },
    ), # with the exception of FGF13: that's not secreted
    af.AnnotDef(
        name = 'gdnf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'GDNF family ligands',
        },
    ), # neurotrophic ligands
    af.AnnotDef(
        name = 'growth_hormone',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Growth hormone family',
        },
    ),
    af.AnnotDef(
        name = 'hedgehog',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hedgehog signaling molecule family',
        },
    ),  # hedgehog proteins are initially membrane bound and can be
        # solubilized later
    af.AnnotDef(
        name = 'hdgf',
        parent = 'ligand',
        resource = {'P51858'},
    ),  # hepatoma-derived growth factor (HDGFL1 is not included because
        # little is known about its role and whether it's secreted)
    af.AnnotDef(
        name = 'igf',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'IGF like family',
        },
    ),
    af.AnnotDef(
        name = 'izumo',
        parent = 'ligand',
        resource = {'Q1ZYL8'},
    ),  # ligands in sperm-egg fusion
    af.AnnotDef(
        name = 'inhibin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Inhibin subunits',
        },
    ),
    af.AnnotDef(
        name = 'interleukin6',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukin 6 type cytokine family',
        },
    ),
    af.AnnotDef(
        name = 'interleukin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interleukins',
        },
    ),
    af.AnnotDef(
        name = 'leucine_rich_glioma_inactivated',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'LGI family',
        },
    ),  # maybe not ligands in a strict sense but don't fit either
        # in other categories
    af.AnnotDef(
        name = 'mia',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'MIA family',
        },
        exclude = {'Q5JRA6', 'Q96PC5'},
    ),
    af.AnnotDef(
        name = 'neuferricin_neudensin',
        parent = 'ligand',
        resource = {'Q8WUJ1', 'Q9UMX5'},
    ),
    af.AnnotDef(
        name = 'netrin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Netrins',
        },
    ),  # secreted axon guidance molecules
    af.AnnotDef(
        name = 'neurotrophin',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurotrophins',
        },
    ),
    af.AnnotDef(
        name = 'oocyte_secreted',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'OOSP family',
        },
    ),  # not sure these are ligands, but this category looks the most likely
        # at least in a recent paper PLAC1 has been described to activate
        # FGFR2 together with FGF7
    af.AnnotDef(
        name = 'prostate_and_testis_expressed',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'PATE family',
        },
    ),  # ligands modulating nicotinic ACh receptors and sperm motility
    af.AnnotDef(
        name = 'pregnancy_specific_glycoprotein',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pregnancy specific glycoproteins',
        },
    ),  # bind to cell surface moieties, regulate other ligands --
        # overall they fit the best to the ligand category
    af.AnnotDef(
        name = 's100_calcium_binding',
        parent = 'ligand',
        resource = {'P31151', 'P80511', 'P05109'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'ligand',
        resource = {'O43866'},
    ),
    af.AnnotDef(
        name = 'tafa_chemokine_like',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'TAFA chemokine like family',
        },
    ),
    af.AnnotDef(
        name = 'prosalusin',
        parent = 'ligand',
        resource = {'Q8N2E6'},
    ),
    af.AnnotDef(
        name = 'transforming_growth_factor_beta',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Transforming growth factor beta family',
        },
    ),
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tumor necrosis factor superfamily',
        },
        exclude = {
            'Q9UNG2', 'P32970', 'O14788', 'P48023',
            'P32971', 'P41273', 'P23510', 'Q06643',
        },
    ),
    af.AnnotDef(
        name = 'vascular_endothelial_growth_factor',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'VEGF family',
        },
    ),
    af.AnnotDef(
        name = 'wnt',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Wnt family',
        },
    ),
    # from Almen 2009
    af.AnnotDef(
        name = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Ligand',
        },
    ),
    af.AnnotDef(
        name = 'delta_like',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Delta',
        },
    ),
    af.AnnotDef(
        name = 'ephrin_b',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'EphB',
        },
    ),
    af.AnnotDef(
        name = 'ig',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'IG_Ligand',
        },
    ),
    af.AnnotDef(
        name = 'jagged',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Jagged',
        },
    ),
    af.AnnotDef(
        name = 'neuroligin',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Neuroligin',
        },
    ),
    af.AnnotDef(
        name = 'nkg2dl',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'NKG2DL',
        },
    ),
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Semaphorins',
        },
    ),
    af.AnnotDef(
        name = 'trem_like',
        parent = 'ligand',
        resource = {'Q6UXN2'},
    ),
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'ligand',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Semaphorin',
        },
    ),
    af.AnnotDef(
        name = 'wnt',
        parent = 'ligand',
        resource = {
            'O96014', 'P56703', 'O00755', 'P04628', 'P56706', 'Q9UBV4',
            'P56705', 'Q9Y6F9', 'O14904', 'P41221', 'Q93098', 'P09544',
            'Q9H1J5', 'Q93097', 'O00744', 'P56704', 'Q9GZT5', 'O14905',
            'Q9H1J7',
        }
    ),
    af.AnnotDef(
        name = 'bmp',
        parent = 'ligand',
        resource = {'P12644', 'P13497'},
    ),
    af.AnnotDef(
        name = 'tgf_beta',
        parent = 'ligand',
        resource = {'P01137'},
    ),
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'ligand',
        resource = {
            'Q14563', 'Q13214', 'Q13275', 'Q99985',
            'Q9NS98', 'O15041', 'O95025',
        }
    ),  # secreted ligands for plexins, regulating axonal growth
    af.AnnotDef(
        name = 'anosmin',
        parent = 'ligand',
        resource = {'P23352'},
    ),

    # ligand regulator
    af.AnnotDef(
        name = 'ligand_regulator',
        scope = 'generic',
        source = 'composite',
        resource = '~ligand_regulator',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'ligand_agonist',
        parent = 'ligand_regulator',
        resource = 'CellChatDB',
        args = {
            'role': 'agonist',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ligand_antagonist',
        parent = 'ligand_regulator',
        resource = 'CellChatDB',
        args = {
            'role': 'antagonist',
        },
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'ligand_regulator',
        resource = {'Q03167'},
    ),
    af.AnnotDef(
        name = 'tgf_beta_binding',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Latent transforming growth factor beta binding proteins'
            ),
        },
    ),
    af.AnnotDef(
        name = 'sparc',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'SPARC family',
        },
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'ligand_regulator',
        resource = {'Q8WTU2'},
    ),
    af.AnnotDef(
        name = 'frizzled_related',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Secreted frizzled-related proteins',
        },
    ),  # secreted proteins binding WNT ligands
    af.AnnotDef(
        name = 'secretoglobin',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Secretoglobins',
        },
    ),  # secreted proteins binding small molecule ligands
    af.AnnotDef(
        name = 'igf_binding',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Insulin like growth factor binding proteins',
        },
    ),
    af.AnnotDef(
        name = 'growth_factor_binder',
        parent = 'ligand_regulator',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Growth factor binding',
        },
    ),
    af.AnnotDef(
        name = 'growth_factor_binder',
        parent = 'ligand_regulator',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Growth Factor-binding',
        },
    ),
    af.AnnotDef(
        name = 'tgf_beta_binding',
        parent = 'ligand_regulator',
        resource = {'Q8N2S1'},
    ),
    af.AnnotDef(
        name = 'glypican',
        parent = 'ligand_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Glypicans',
        },
    ),
    af.AnnotDef(
        name = 'glypican',
        parent = 'ligand_regulator',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Glypican',
        },
    ),
    af.AnnotDef(
        name = 'syndecan',
        parent = 'ligand_regulator',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Syndecan',
        },
    ),
    # growth factor binder or regulator
    af.AnnotDef(
        name = 'growth_factor_binder',
        parent = 'ligand_regulator',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'growth factor binding',
        },
    ),
    af.AnnotDef(
        name = 'ligand_regulator',
        scope = 'generic',
        resource = 'Matrisome',
        args = {
            'mainclass': 'Matrisome-associated',
            'subclass': 'Secreted Factors',
        },
        limit = ('growth_factor_binder', 'secreted_enzyme'),
        enabled = False,
    ),  # to be checked later


    # cell surface ligand
    af.AnnotDef(
        name = 'cell_surface_ligand',
        scope = 'generic',
        source = 'composite',
        resource = '~cell_surface_ligand',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'cell_surface_ligand',
        scope = 'generic',
        resource = 'connectomeDB2020',
        args = {
            'role': 'ligand',
            'location': 'plasma membrane',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_ligand',
        resource = 'CellChatDB',
        scope = 'generic',
        args = {
            'role': 'ligand',
            'category': 'Cell-Cell Contact',
        },
    ),
    af.AnnotDef(
        name = 'cell_surface_ligand',
        resource = 'Baccin2019',
        scope = 'generic',
        args = {
            'mainclass': 'ligand',
            'location': {'membrane', 'both'},
        },
        # `exclude` defined in `excludes` dict
    ),  # quite some wrong annotations, also many receptors
    af.AnnotDef(
        name = 'tumor_necrosis_factor',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tumor necrosis factor superfamily',
        },
        exclude = {'O75888'},
    ),
    af.AnnotDef(
        name = 'cell_surface_ligand',
        scope = 'generic',
        resource_name = 'GO_Intercell',
        resource = af.AnnotOp(
            annots = (
                'cell_surface',
                'ligand_go',
            ),
            op = set.intersection,
        ),
        enabled = False, # to be checked, disabled until then
    ),
    af.AnnotDef(
        name = 'cell_surface_ligand',
        scope = 'generic',
        resource = 'CellPhoneDB',
        args = {
            'method': lambda a: (
                not a.receptor and (
                    a.peripheral or
                    a.transmembrane
                )
            ),
        },
    ),
    # surface ligand subclasses from HGNC
    af.AnnotDef(
        name = 'b7_family',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'B7 family',
        },
    ),
    af.AnnotDef(
        name = 'butyrophilin',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Butyrophilins',
        },
    ),
    af.AnnotDef(
        name = 'ephrin',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ephrins',
        },
    ),
    af.AnnotDef(
        name = 'neuregulin',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuregulins',
        },
    ),  # ligands for various ERBB receptors
    af.AnnotDef(
        name = 'hedgehog',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hedgehog signaling molecule family',
        },
    ),  # hedgehog proteins are initially membrane bound and can be
        # solubilized later
    af.AnnotDef(
        name = 'izumo',
        parent = 'cell_surface_ligand',
        resource = {'Q8IYV9', 'Q6UXV1', 'Q5VZ72'},
    ),  # ligands in sperm-egg fusion
    af.AnnotDef(
        name = 'nectin',
        parent = 'cell_surface_ligand',
        resource = {'O95727', 'Q15223'},
    ),  # ligands for T-lymphocytes
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Semaphorins',
        },
        exclude = {
            'Q14563', 'Q13214', 'Q13275', 'Q99985',
            'Q9NS98', 'O15041', 'O95025',
        },
    ),  # surface bound ligands for plexins, regulating axonal growth
    af.AnnotDef(
        name = 'anosmin',
        parent = 'cell_surface_ligand',
        resource = {'P23352'},
    ),
    af.AnnotDef(
        name = 'mhc',
        parent = 'cell_surface_ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'MHC',
        },
    ),
    af.AnnotDef(
        name = 'semaphorin',
        parent = 'cell_surface_ligand',
        resource = 'Almen2009',
        args = {
            'classes': 'Semaphorins',
        },
    ),  # surface bound ligands for plexins, regulating axonal growth


    # adhesion (generic, cell-cell or cell-matrix undistinguished)
    af.AnnotDef(
        name = 'adhesion',
        scope = 'generic',
        source = 'composite',
        resource = af.AnnotOp(
            annots = (
                '~adhesion',
                '~cell_adhesion',
                '~matrix_adhesion',
            ),
            op = set.union,
        ),
        transmitter = True,
        receiver = True,
        # transmitter & receiver are not valid for all members of this
        # category but we have no better option now than keep all these
        # proteins as potential transmitters and receivers
    ),
    af.AnnotDef(
        name = 'adhesion',
        resource = af.AnnotOp(
            annots = (
                '~cell_adhesion~HGNC',
                '~matrix_adhesion~HGNC',
            ),
            op = set.union,
        ),
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'adhesion',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Cell adhesion',
        },
        limit = 'cell_surface',
    ),  # with limiting to the cell surface, it's a nice
        # collecion of adhesion proteins (267)
    af.AnnotDef(
        name = 'adhesion',
        scope = 'generic',
        resource = 'CellPhoneDB',
        args = {
            'integrin': bool,
        },
    ),
    af.AnnotDef(
        name = 'adhesion',
        scope = 'generic',
        resource = 'MCAM',
    ),
    af.AnnotDef(
        name = 'adhesion',
        scope = 'generic',
        resource = 'Adhesome',
        args = {
            'mainclass': 'Adhesion receptor',
        },
    ),  # both cell-cell and cell-matrix adhesion, good as it is
    af.AnnotDef(
        name = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Adhesion',
        },
    ),
    af.AnnotDef(
        name = 'ig_like',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'IG_AdhesionProteins',
        },
    ),
    af.AnnotDef(
        name = 'mpz',
        parent = 'adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'IG_MPZ',
        },
    ),
    af.AnnotDef(
        name = 'ly6_plaur',
        parent = 'adhesion',
        resource = {'O95274', 'Q8N6Q3', 'Q8TDM5', 'Q9BY14', 'Q17RY6'},
    ),  # both cell and matrix adhesion

    # cell-cell adhesion
    af.AnnotDef(
        name = 'cell_adhesion',
        scope = 'generic',
        source = 'composite',
        resource = af.AnnotOp(
            annots = '~cell_adhesion',
            op = set.union,
        ),
        transmitter = True,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'cell_adhesion',
        scope = 'generic',
        resource = 'Cellinker',
        args = {
            'type': 'Cell adhesion',
        },
        exclude = {
            'Q02297', 'Q14952', 'P29965', 'Q14943', 'Q5JRA6', 'P35221',
            'P35222', 'Q9UI47', 'Q9UQB3', 'Q9UBT7', 'Q14953', 'Q5T4B2',
            'Q8N743', 'Q13895', 'Q7L5Y9', 'Q99689', 'Q8NHK3',
        },
    ),
    af.AnnotDef(
        name = 'cell_adhesion',
        scope = 'generic',
        resource = 'Zhong2015',
        args = {
            'type': {
                'cell-cell adhesion',
                'iCAM',
                'myelin interactions',
            },
        },
        exclude = {
            'Q02297', 'Q14952', 'P29965', 'Q14943', 'Q5JRA6', 'P35221',
            'P35222', 'Q9UI47', 'Q9UQB3', 'Q9UBT7', 'Q14953', 'Q5T4B2',
            'Q8N743', 'Q13895', 'Q7L5Y9', 'Q99689', 'Q8NHK3',
        },
    ),
    af.AnnotDef(
        name = 'icam',
        parent = 'cell_adhesion',
        resource = 'Zhong2015',
        args = {
            'type': 'iCAM',
        },
    ),
    af.AnnotDef(
        name = 'myelin_adhesion',
        parent = 'cell_adhesion',
        resource = 'Zhong2015',
        args = {
            'type': 'myelin interactions',
        },
    ),
    af.AnnotDef(
        name = 'cell_adhesion',
        resource = 'GO_Intercell',
        scope = 'generic',
        args = {
            'mainclass': 'adhesion to other cells',
        },
        limit = 'cell_surface',
    ),
    af.AnnotDef(
        name = 'cell_adhesion',
        resource = '~cell_adhesion~Almen2009',
        resource_name = 'Almen2009',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'beta_protocadherin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'ProtocadherinsBeta',
        },
    ),  # cell-cell adhesion especially between neurons
    af.AnnotDef(
        name = 'protocadherin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'ProtocadherinsOther',
        },
    ),  # cell-cell adhesion
    af.AnnotDef(
        name = 'classical_cadherin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'CadherinClassic',
        },
    ),  # cell-cell adhesion
    af.AnnotDef(
        name = 'desmosomal_cadherin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'CadherinOther',
        },
    ),  # cell-cell adhesion
    # cell-cell adhesion HGNC
    af.AnnotDef(
        name = 'cell_adhesion',
        resource = '~cell_adhesion~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = '7d_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': '7D cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'ceacam',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Carcinoembryonic antigen related '
                'cell adhesion molecule family'
            ),
        },
    ), # in plasma membrane; cell-cell adhesion and receptors
    af.AnnotDef(
        name = 'cadherin_related',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cadherin related',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'major_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Major cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'non_clustered_protocadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Non-clustered protocadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'type1_classical_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type I classical cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'type2_classical_cadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Type II classical cadherins',
        },
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'nectin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Nectins and nectin-like molecules',
        },
        exclude = {'O95727', 'Q15223'},
    ), # cell-cell adhesion
    af.AnnotDef(
        name = 'neurexin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neurexins',
        },
    ), # cell-cell adhesion for neurons
    af.AnnotDef(
        name = 'neurexin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Neurexin',
        },
    ),  # these are also receptors
    af.AnnotDef(
        name = 'neuroligin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Neuroligins',
        },
    ), # cell-cell adhesion for neurons
    af.AnnotDef(
        name = 'neuroligin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Neuroligin',
        },
    ), # cell-cell adhesion for neurons; these are also ligands
    af.AnnotDef(
        name = 'clarin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Clarins',
        },
    ),  # regulation of cell-cell adhesion and synapsis in ear and retina
    af.AnnotDef(
        name = 'protocadherin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Clustered protocadherins',
        },
    ),  # cell-cell adhesion in brain neuronal connections
    af.AnnotDef(
        name = 'ig_like',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ig-like cell adhesion molecule family',
        },
        exclude = {'Q9HCN6', 'Q14CZ8'},
    ),
    af.AnnotDef(
        name = 'igcam_cxadr_like',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'IgCAM CXADR-related subfamily',
        },
    ),
    af.AnnotDef(
        name = 'iglon',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'IgLON cell adhesion molecules',
        },
        exclude = {'A6NGN9'},
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'cell_adhesion',
        resource = {
            'P23229', 'Q13349', 'Q13797', 'P20701', 'P38570', 'P05107',
            'P26010',
        },
    ),
    af.AnnotDef(
        name = 'receptor_tyrosine_phosphatase',
        parent = 'cell_adhesion',
        resource = {'P28827', 'O14522'},
    ),  # most of the PTPRs are not adhesion molecules but only adhesion
        # receptors or other receptors
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'cell_adhesion',
        resource = {'P16070'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'cell_adhesion',
        resource = {'P30203'},
    ),
    af.AnnotDef(
        name = 'selectin',
        parent = 'cell_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Selectins',
        },
    ),
    af.AnnotDef(
        name = 'selectin',
        parent = 'cell_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Selectin',
        },
    ),  # cell-cell adhesion between immune cells

    # cell-matrix adhesion
    af.AnnotDef(
        name = 'matrix_adhesion',
        scope = 'generic',
        source = 'composite',
        resource = af.AnnotOp(
            annots = '~matrix_adhesion',
            op = set.union,
        ),
        transmitter = False,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'matrix_adhesion',
        scope = 'generic',
        resource = 'Cellinker',
        args = {
            'role': 'receptor',
            'type': {
                'ECM-receptor interaction',
            },
        },
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'matrix_adhesion',
        resource = 'Integrins',
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'matrix_adhesion',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Integrin',
        },
    ),
    af.AnnotDef(
        name = 'matrix_adhesion',
        scope = 'generic',
        resource = 'Zhong2015',
        args = {
            'type': {
                'matrix adhesion',
                'focal adhesion',
            },
        },
        exclude = {'Q14511', 'O60711', 'Q7Z4I7', 'Q96AC1', 'Q9H792'},
    ),
    af.AnnotDef(
        name = 'focal_adhesion',
        parent = 'matrix_adhesion',
        resource = 'Zhong2015',
        args = {
            'type': 'focal adhesion',
        },
    ),
    af.AnnotDef(
        name = 'matrix_adhesion',
        resource = 'GO_Intercell',
        scope = 'generic',
        args = {
            'mainclass': 'adhesion to matrix',
        },
        limit = 'cell_surface',
    ),
    af.AnnotDef(
        name = 'focal_adhesion',
        parent = 'matrix_adhesion',
        resource = 'Ramilowski_location',
        args = {
            'location': 'focal adhesion',
        },
        limit = 'cell_surface',
        exclude = {
            'O14976', 'O15231', 'O43294', 'O60711', 'P07332', 'P50552',
            'P56945', 'Q13153', 'Q14451', 'Q14511', 'Q155Q3', 'Q15654',
            'Q7Z4I7', 'Q86TP1', 'Q8IVT2', 'Q8IZW8', 'Q8N264', 'Q8WX93',
            'Q92502', 'Q96QB1', 'Q9H792', 'Q9HBI0', 'Q9HBI1', 'Q9NQ75',
            'Q9UGI8', 'Q9UGP4',
        },
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'matrix_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Integrin',
        },
    ),  # matrix adhesion; these are also receptors
    af.AnnotDef(
        name = 'sarcoglycan',
        parent = 'matrix_adhesion',
        resource = 'Almen2009',
        args = {
            'classes': 'Sarcoglycan',
        },
    ),  # cell-matrix adhesion for muscle cells
    af.AnnotDef(
        name = 'matrix_adhesion',
        resource = '~matrix_adhesion~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'adhesion_gprotein_coupled_receptor',
        parent = 'matrix_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Adhesion G protein-coupled receptors, subfamily A',
                'Adhesion G protein-coupled receptors, subfamily B',
                'Adhesion G protein-coupled receptors, subfamily C',
                'Adhesion G protein-coupled receptors, subfamily D',
                'Adhesion G protein-coupled receptors, subfamily E',
                'Adhesion G protein-coupled receptors, subfamily F',
                'Adhesion G protein-coupled receptors, subfamily G',
                'Adhesion G protein-coupled receptors, subfamily L',
                'Adhesion G protein-coupled receptors, subfamily V',
            },
        }, # cell-matrix adhesion
    ),
    af.AnnotDef(
        name = 'matrix_adhesion',
        parent = 'adhesion',
        resource = {'Q9HCN6', 'Q14CZ8'},
    ),
    af.AnnotDef(
        name = 'integrin',
        parent = 'matrix_adhesion',
        resource = {
            'P06756', 'Q9UKX5', 'P08648', 'P11215', 'P26006', 'Q13683',
            'P20702', 'O75578', 'P13612', 'P17301', 'P56199', 'P08514',
            'P18564', 'O95965', 'P18084', 'P05556', 'P26012', 'P16144',
            'P05106',
        },
    ),
    af.AnnotDef(
        name = 'mucin',
        parent = 'matrix_adhesion',
        resource = af.AnnotOp(
            annots = (
                af.AnnotDef(
                    name = 'mucin_hgnc',
                    resource ='HGNC',
                    args = {
                        'mainclass': 'Mucins',
                    },
                ),
                'transmembrane'
            ),
            op = set.intersection
        ),
    ),  # membrane bound mucins
    af.AnnotDef(
        name = 'proteoglycan',
        parent = 'matrix_adhesion',
        resource = {'Q6UVK1'},
    ),

    # cell-matrix adhesion regulators
    af.AnnotDef(
        name = 'matrix_adhesion_regulator',
        resource = '~matrix_adhesion_regulator',
        scope = 'generic',
        source = 'composite',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'matrix_adhesion_regulator',
        resource = {'Q08380'},
    ),

    # surface enzyme
    af.AnnotDef(
        name = 'cell_surface_enzyme',
        scope = 'generic',
        source = 'composite',
        resource = af.AnnotOp(
            annots = (
                '~cell_surface_enzyme',
                '~cell_surface_peptidase',
            ),
            op = set.union,
        ),
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'cell_surface_enzyme',
        resource = 'GO_Intercell',
        scope = 'generic',
        args = {
            'mainclass': 'enzyme',
        },
        limit = 'cell_surface',
        avoid = 'receptor',
    ),  # to be checked
    af.AnnotDef(
        name = 'cell_surface_enzyme',
        resource = 'Surfaceome',
        scope = 'generic',
        args = {
            'mainclass': 'Enzymes',
        },
        exclude = {
            'Q8TCJ2', 'Q9NPH5', 'Q9UKF2', 'P04843', 'P19021', 'Q9UK23',
            'P52961', 'P09958', 'P11117', 'Q96JJ7', 'P40126', 'P08842',
            'Q13444', 'Q8WUD6', 'O00391', 'Q15125',
        },
    ),  # looks all right
    af.AnnotDef(
        name = 'cell_surface_enzyme',
        resource = af.AnnotOp(
            annots = (
                '~cell_surface_enzyme~HGNC',
                '~cell_surface_peptidase~HGNC',
            ),
            op = set.union,
        ),
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'ectonucleotide_phosphatase',
        parent = 'cell_surface_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Ectonucleotide pyrophosphatase/phosphodiesterase family'
            ),
        },
    ),  # maybe not all bound to the surface but most of them
    af.AnnotDef(
        name = 'hyaluronidase',
        parent = 'cell_surface_enzyme',
        resource = {'Q9UHN6', 'Q12891', 'P38567', 'Q2M3T9'},
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'cell_surface_enzyme',
        resource = {
            'P98073', 'Q9Y5Q5', 'Q9BYE2', 'Q9H3S3', 'O15393', 'P05981',
            'Q9NRS4',
        },
    ),

    # cell surface peptidase (protease)
    af.AnnotDef(
        name = 'cell_surface_peptidase',
        scope = 'generic',
        source = 'composite',
        resource = '~cell_surface_peptidase',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'cell_surface_peptidase',
        resource = '~cell_surface_peptidase~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'insulin_degrading',
        parent = 'cell_surface_peptidase',
        resource = {'P14735'},
    ),
    af.AnnotDef(
        name = 'm1_metallopeptidase',
        parent = 'cell_surface_peptidase',
        resource = {'Q6Q4G3', 'Q9UIQ6', 'Q07075', 'P15144', 'Q9UKU6'},
    ),  # cleave mostly peptide ligands, hormones like TRH, angiotensin, etc
    af.AnnotDef(
        name = 'm16_metallopeptidase',
        parent = 'cell_surface_peptidase',
        resource = {'P14735'},
    ),  # acts on peptide hormones
    af.AnnotDef(
        name = 'm10_metallopeptidase',
        parent = 'cell_surface_peptidase',
        resource = {'P51511', 'P51512', 'Q9ULZ9', 'Q9Y5R2', 'Q9NPA2'},
    ),
    af.AnnotDef(
        name = 'm13_metallopeptidase',
        parent = 'cell_surface_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'M13 metallopeptidases',
        },
    ),
    af.AnnotDef(
        name = 'm14_carboxypeptidase',
        parent = 'cell_surface_peptidase',
        resource = {'P14384', 'O75976', 'Q8IVL8', },
    ),
    af.AnnotDef(
        name = 'vanin',
        parent = 'cell_surface_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Vanins',
        },
        exclude = {'P43251'},
    ),
    af.AnnotDef(
        name = 'serine_protease',
        parent = 'cell_surface_peptidase',
        resource = {
            'Q8IU80', 'Q9Y5Q5', 'Q7RTY9', 'Q9BYE2', 'Q6ZMR5', 'Q9H3S3',
            'Q86WS5', 'Q9Y6M0', 'Q9Y5Y6', 'O60235', 'Q9NRR2', 'Q7Z410',
            'Q6ZWK6', 'Q7RTY8', 'P05981', 'Q86T26', 'Q6UWB4', 'Q9NRS4',
            'A4D1T9', 'O15393', 'Q92743', 'Q9UL52', 'Q16651'
        },
    ),


    # secreted enzyme
    af.AnnotDef(
        name = 'secreted_enzyme',
        scope = 'generic',
        source = 'composite',
        resource = af.AnnotOp(
            annots = (
                '~secreted_enzyme',
                '~secreted_peptidase',
            ),
            op = set.union,
        ),
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'secreted_enzyme',
        resource = 'GO_Intercell',
        scope = 'generic',
        args = {
            'mainclass': 'enzyme',
        },
        limit = 'secreted',
        avoid = ('receptor', 'ligand', 'plasma_membrane_transmembrane'),
    ),
    af.AnnotDef(
        name = 'secreted_enzyme',
        resource = af.AnnotOp(
            annots = (
                '~secreted_enzyme~HGNC',
                '~secreted_peptidase~HGNC',
            ),
            op = set.union,
        ),
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'phospholipase',
        parent = 'secreted_enzyme',
        resource = {
            'Q53H76', 'Q8NCC3', 'Q9NZK7', 'Q9BX93', 'P04054', 'Q13093',
            'Q5R387', 'Q9NZ20', 'P14555', 'Q9BZM1', 'Q9UNK4', 'O15496',
            'Q9BZM2', 'P39877',
        },
    ),  # secreted enzymes acting on phospholipids
    af.AnnotDef(
        name = 'defensin',
        parent = 'secreted_enyzme',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Defensins, alpha',
                'Defensins, beta',
            },
        },
    ),  # permeabilizing microorganism membranes or
        # binding to microorganism surfaces
    af.AnnotDef(
        name = 'lysozym',
        parent = 'secreted_enyzme',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Lysozymes, c-type',
                'Lysozymes, g-type'
            },
        },
    ),  # bacteriolytic proteins, some involved in sperm-egg fertilization
    af.AnnotDef(
        name = 'galactosidase',
        parent = 'secreted_enyzme',
        resource = {'Q6UWU2', 'Q8IW92'},
    ),  # secreted galactosidases
    af.AnnotDef(
        name = 'lipase',
        parent = 'secreted_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Lipases',
        },
        limit = 'secreted',
    ),  # secreted lipases
    af.AnnotDef(
        name = 'paraoxonase',
        parent = 'secreted_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Paraoxonases',
        },
    ),  # secreted enzymes hydrolysing lactons and other metabolites
    af.AnnotDef(
        name = 'lipocalin',
        parent = 'secreted_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Lipocalins',
        },
        limit = 'secreted',
    ),  # secreted lipases
    af.AnnotDef(
        name = 'lipase',
        parent = 'secreted_enzyme',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Lipid degradation',
        },
        limit = 'secreted',
    ),
    af.AnnotDef(
        name = 'scavenger_receptor_cysteine_rich',
        parent = 'secreted_enzyme',
        resource = {
            'P58215', 'Q96JB6', 'P05156', 'P56730', 'Q96JK4', 'Q9Y4K0',
        },
    ),
    af.AnnotDef(
        name = 'glutathione_peroxidase',
        parent = 'secreted_enzyme',
        resource = {'P59796', 'P49908', 'P22352'},
    ),
    af.AnnotDef(
        name = 'biotinidase',
        parent = 'secreted_enzyme',
        resource = {'P43251'},
    ),
    af.AnnotDef(
        name = 'hyaluronidase',
        parent = 'secreted_enzyme',
        resource = {'Q12794', 'Q8WUJ3', 'O43820'},
    ),
    af.AnnotDef(
        name = 'ribonuclease',
        parent = 'secreted_enzyme',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ribonuclease A family',
        },
        exclude = {'P10153', 'P03950'},
    ),

    # secreted peptidase
    af.AnnotDef(
        name = 'secreted_peptidase',
        resource = '~secreted_peptidase',
        source = 'composite',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'secreted_peptidase',
        scope = 'generic',
        resource = 'GO_Intercell',
        limit = 'secreted',
        args = {
            'mainclass': 'peptidase',
        },
    ),
    af.AnnotDef(
        name = 'secreted_peptidase',
        resource = 'UniProt_keyword',
        scope = 'generic',
        args = {
            'keyword': {
                'Protease',
                'Serine protease',
            },
        },
        limit = 'secreted',
    ),  # looks all right
    af.AnnotDef(
        name = 'collagen_degrading',
        parent = 'secreted_peptidase',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Collagen degradation',
        },
    ),  # very good
    af.AnnotDef(
        name = 'cathepsin',
        parent = 'secreted_peptidase',
        resource = {'P08311', 'P07711', 'P43235', 'P25774', 'P07339'},
    ),
    af.AnnotDef(
        name = 'secreted_peptidase',
        resource = '~secreted_peptidase~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    # subclasses from HGNC
    af.AnnotDef(
        name = 'adamts',
        parent = 'secreted_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'ADAM metallopeptidases with thrombospondin type 1 motif'
            ),
        },
    ),
    af.AnnotDef(
        name = 'adamts_like',
        parent = 'secreted_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'ADAMTS like',
        },
    ),
    af.AnnotDef(
        name = 'heparanase',
        parent = 'secreted_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Heparanases',
        },
    ),  # act on heparin and heparane-sulphate
    af.AnnotDef(
        name = 'm14_carboxypeptidase',
        parent = 'secreted_peptidase',
        resource = {
            'P16870', 'P15169', 'Q8IUX7', 'Q66K79', 'P48052', 'Q96SM3',
            'Q8WXQ8', 'P15086', 'Q96IY4', 'P15085', 'Q8N4T0', 'Q9HB40',
            'P22792', 'Q9UI42', 'Q8N436', 'Q9Y646',
        },
    ),
    af.AnnotDef(
        name = 'm1_metallopeptidase',
        parent = 'secreted_peptidase',
        resource = {'Q9H4A4'},
    ),  #
    af.AnnotDef(
        name = 'm16_metallopeptidase',
        parent = 'secreted_peptidase',
        resource = {'P14735'},
    ),  # acts on peptide hormones
    af.AnnotDef(
        name = 'immune_serin_protease',
        parent = 'secreted_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Granule associated serine proteases of immune defence'
            ),
        },
    ),  # secreted by granulocytes as part of the immune response
    af.AnnotDef(
        name = 'm10_metallopeptidase',
        parent = 'secreted_peptidase',
        resource = {
            'Q9H239', 'P09237', 'P09238', 'P03956', 'P08253', 'P24347',
            'P39900', 'P45452', 'Q9NRE1', 'P22894', 'Q99542', 'Q8N119',
            'O60882', 'P14780', 'P08254',
        },
    ),  # secreted matrix metallopeptidases, many act on the ECM
    # specific classes from HGNC
    af.AnnotDef(
        name = 'chymotrypsin_like_elastase',
        parent = 'extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chymotrypsin like elastases',
        },
    ),  # involved in ECM dynamics and remodeling
    af.AnnotDef(
        name = 'kallikrein',
        parent = 'extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Kallikreins',
        },
    ),  # extracellular serine proteases, involved in ECM dynamics
    af.AnnotDef(
        name = 'pappalysin',
        parent = 'extracellular_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pappalysins',
        },
    ),  # cleave IGFBPs
    af.AnnotDef(
        name = 'serine_protease',
        parent = 'secreted_peptidase',
        resource = 'HGNC',
        args = {
            'mainclass': 'Serine proteases',
        },
        exclude = {
            'Q9UHE8', 'P36776', 'Q9NQE7', 'Q6UWY2', 'Q7Z5A4', 'P57727',
            'Q8IU80', 'Q9Y5Q5', 'Q7RTY9', 'Q9BYE2', 'Q6ZMR5', 'Q9H3S3',
            'Q86WS5', 'Q9Y6M0', 'Q9Y5Y6', 'O60235', 'Q9NRR2', 'Q7Z410',
            'Q6ZWK6', 'Q7RTY8', 'P05981', 'Q86T26', 'Q6UWB4', 'Q9NRS4',
            'O43464', 'Q9UI38',
        },
    ),
    af.AnnotDef(
        name = 'matrix_matalloproteinase',
        parent = 'secreted_peptidase',
        resource = {
            'Q9NRE1', 'Q9NPA2', 'P14780', 'Q9H239', 'P22894', 'P08253',
            'P09237', 'P39900', 'P45452', 'Q9ULZ9', 'P51512', 'P03956',
        }
    ),
    af.AnnotDef(
        name = 'adamts',
        parent = 'secreted_peptidase',
        resource = {'Q8N6G6', 'Q6UY14'},
    ),
    af.AnnotDef(
        name = 'stromelysin',
        parent = 'secreted_peptidase',
        resource = {'P08254', 'P09238'},
    ),
    af.AnnotDef(
        name = 'carboxypeptidase',
        parent = 'secreted_peptidase',
        resource = {'Q66K79'},
    ),
    af.AnnotDef(
        name = 'polyserase',
        parent = 'secreted_peptidase',
        resource = {'Q5K4E3'},
    ),


    # secreted peptidase inhibitor
    af.AnnotDef(
        name = 'secreted_peptidase_inhibitor',
        scope = 'generic',
        source = 'composite',
        resource = '~secreted_peptidase_inhibitor',
    ),
    af.AnnotDef(
        name = 'serpin',
        parent = 'secreted_peptidase_inhibitor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Serpin peptidase inhibitors',
        },
        limit = 'secreted',
    ),
    af.AnnotDef(
        name = 'inter_alpha_trypsin_inhibitor',
        parent = 'secreted_peptidase_inhibitor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Inter-alpha-trypsin inhibitor heavy chains',
        },
    ),  # protease inhibitors in plasma
    af.AnnotDef(
        name = 'secreted_peptidase_inhibitor',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': {
                'Protease inhibitor',
                'Serine protease inhibitor',
            },
        },
        limit = 'secreted',
    ),  # looks all right
    af.AnnotDef(
        name = 'cystatin',
        parent = 'secreted_peptidase_inhibitor',
        resource = 'Matrisome',
        args = {
            'subsubclass': 'Cystatin',
        },
        exclude = {'P04080', 'P01040'},
    ),
    af.AnnotDef(
        name = 'serpin',
        parent = 'secreted_peptidase_inhibitor',
        resource = {'P01009'},
    ),
    af.AnnotDef(
        name = 'tissue_inhibitor_of_metallopeptidases',
        parent = 'secreted_peptidase_inhibitor',
        resource = {'P35625'},
    ),
    af.AnnotDef(
        name = 'serine_peptidase_inhibitor',
        parent = 'secreted_peptidase_inhibitor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Serine peptidase inhibitors, Kazal type',
        },
    ),
    af.AnnotDef(
        name = 'wap_four_disulfide_core_domain_containing',
        parent = 'secreted_peptidase_inhibitor',
        resource = 'HGNC',
        args = {
            'mainclass': 'WAP four-disulfide core domain containing',
        },
        exclude = {'P23352'},
    ),
    af.AnnotDef(
        name = 'tissue_inhibitor_of_metallopeptidases',
        parent = 'secreted_peptidase_inhibitor',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tissue inhibitor of metallopeptidases',
        },
    ),


    # transporter
    af.AnnotDef(
        name = 'transporter',
        resource = af.AnnotOp(
            annots = (
                '~transporter',
                '~ion_channel',
            ),
            op = set.union,
        ),
        scope = 'generic',
        source = 'composite',
        transmitter = False,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'transporter',
        resource = 'Surfaceome',
        scope = 'generic',
        args = {
            'mainclass': 'Transporters',
        },
        exclude = {
            'Q7L1I2', 'Q9ULQ1', 'Q05940', 'Q9BZC7', 'Q8NBW4', 'P54219',
            'Q9P2U8', 'Q8IY34', 'Q8TED4', 'Q9UN42', 'Q9P2U7', 'Q8NCC5',
            'Q9H598', 'Q8NHS3', 'Q9NRX5', 'Q9H1V8', 'Q496J9', 'Q6J4K2',
            'Q96T83', 'Q9NP78', 'A6NFC5', 'Q8TBB6', 'O00400', 'Q8WWZ7',
            'Q71RS6', 'Q9GZU1', 'O95528', 'Q8NDX2', 'O43826', 'O94778',
            'Q9HD20', 'Q9UGQ3',
        },
    ),  # some intracellular transporters added to exclude
    af.AnnotDef(
        name = 'transporter',
        scope = 'generic',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'transport',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'transporter',
        resource = '~transporter~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'transporter',
        parent = 'transporter',
        resource_name = 'DGIdb',
        resource = '~transporter~DGIdb',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'abc',
        parent = 'transporter',
        resource = 'DGIdb',
        args = {
            'category': 'ABC TRANSPORTER',
        },
        exclude = {
            'P28288', 'Q93050', 'P33897', 'P54652', 'P21283', 'Q06055',
            'Q9BRX2', 'Q9NP78', 'P24539', 'P24539', 'Q92736', 'Q9UBJ2',
            'O14678', 'P78363', 'Q9NUT2', 'P06576', 'P62328', 'O14983',
            'Q96LB4', 'Q99437', 'Q9NRK6', 'O75964', 'P48047', 'Q5VTU8',
            'P98194', 'P21281', 'Q13488', 'Q9NR96', 'Q52LC2', 'Q15904',
            'P38606', 'Q9BZC7', 'Q7Z4Y8', 'P26678', 'P36542', 'P48201',
            'P16615', 'P30405', 'O15533', 'Q8WWZ7', 'P50570', 'O75027',
            'O75947', 'O75110', 'P25705', 'P56381', 'Q93084', 'Q03519',
            'O94823', 'P21917', 'Q8N8Y2', 'Q03518', 'P35462', 'P61421',
            'O60423', 'O43861', 'P35670', 'O00631', 'P30049', 'O75534',
            'Q16864', 'P63165', 'Q8NHE4', 'Q9ULM6', 'Q9UI12', 'Q96A05',
            'Q9HD20', 'Q9UHG3', 'P05496', 'P27449', 'Q9Y2G3',
        },
    ),
    af.AnnotDef(
        name = 'transporter',
        scope = 'generic',
        resource = 'Almen2009',
        parent = 'transporter',
        args = {
            'mainclass': 'Transporters',
        },
        limit = 'plasma_membrane_transmembrane',
        exclude = {
            'P32856', 'O95183', 'Q13277', 'Q9UNK0',
        }
    ),
    af.AnnotDef(
        name = 'transporter',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Transport',
        },
        limit = 'plasma_membrane_transmembrane',
    ),

    # transporters from HGNC
    af.AnnotDef(
        name = 'aquaporin',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'Aquaporins',
        },
    ),
    af.AnnotDef(
        name = 'bestrophin',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'Bestrophins',
        },
    ),
    af.AnnotDef(
        name = 'abcc',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'ATP binding cassette subfamily C',
        },
    ),
    af.AnnotDef(
        name = 'abcg',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'ATP binding cassette subfamily G',
        },
    ),
    af.AnnotDef(
        name = 'hk_atpase',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'ATPase H+/K+ transporting',
        },
    ),
    af.AnnotDef(
        name = 'sodium_potassium_atpase',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'ATPase Na+/K+ transporting subunits',
                'Na+/K+ transporting ATPase interacting',
            },
        },
    ),
    af.AnnotDef(
        name = 'cnnm_metal',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'Cyclin and CBS domain divalent metal cation '
                'transport mediators'
            ),
        },
    ),
    af.AnnotDef(
        name = 'pannexin',
        parent = 'transporter',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pannexins',
        },
    ),  # as half channels they release ATP, Ca and other substances to
        # the extracellular space

    # specific subclasses from Almen 2009
    af.AnnotDef(
        name = 'abca',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'ABCA',
        },
    ),  # lipid transporters
    af.AnnotDef(
        name = 'abcc',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'ABCC',
        },
    ),
    af.AnnotDef(
        name = 'abcg',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'ABCG',
        },
    ),
    af.AnnotDef(
        name = 'solute_carrier',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': {
                'AMAC', 'APC', 'SLC1', 'SLC10',
                'SLC11', 'SLC12', 'SLC13', 'SLC14',
                'SLC15', 'SLC16', 'SLC17', 'SLC18',
                'SLC19', 'SLC2', 'SLC20', 'SLC22',
                'SLC23', 'SLC24', 'SLC26', 'SLC27',
                'SLC28', 'SLC29', 'SLC3', 'SLC30',
                'SLC31', 'SLC34', 'SLC36', 'SLC38',
                'SLC39', 'SLC4', 'SLC40', 'SLC41',
                'SLC42(Rh)', 'SLC43', 'SLC44', 'SLC45',
                'SLC46', 'SLC5', 'SLC6', 'SLC7',
                'SLC8', 'SLC9', 'SLCO',
            },
        },
        exclude = {
            'Q8TBB6', 'Q8TE54', 'Q8IY34', 'Q9P2U8', 'Q9P2U7', 'Q8NDX2',
            'Q8NHS3', 'Q05940', 'P54219', 'Q9UGQ3', 'O95528', 'Q8N4V2',
            'Q6J4K2', 'Q71RS6', 'Q8TE54', 'Q6P1M0', 'Q6PML9', 'Q8TAD4',
            'Q6NXT4', 'O14863', 'Q8NEW0', 'Q99726', 'Q9BRI3', 'Q8NBW4',
            'Q92504', 'Q9C0K1', 'Q96H72', 'Q9UMX9', 'Q9H1V8', 'Q8TBB6',
            'Q8IVB4', 'Q9Y2E8',
        },
    ),  # transporters for various compounds, e.g. amino acids, bile acids,
        # metal ions, other inorganic and organic ions, urea, oligopeptides,
        # vitamins, sugars, steroids, organic acids, fatty acids,
        # pyrimidines, purines, nucleosides, 
        # we can split this group later
    af.AnnotDef(
        name = 'magnesium',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'NIPA',
        },
    ),
    af.AnnotDef(
        name = 'sodium_potassium_atpase',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'NKAIN',
        },
    ),
    af.AnnotDef(
        name = 'sphingolipid',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'Spinster',
        },
    ),
    af.AnnotDef(
        name = 'xk',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'XK',
        },
    ),  # transporters for amino acids,
        # scramblases for phosphatidylserine,
        # and who knows what else
    af.AnnotDef(
        name = 'aquaporin',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'Aquaporins',
        },
    ),
    af.AnnotDef(
        name = 'auxiliary_transport_unit',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'AuxillaryTransportUnit',
        },
        exclude = {'Q9UN42'},
    ),
    af.AnnotDef(
        name = 'bestrophin',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'Bestrophin',
        },
    ),
    af.AnnotDef(
        name = 'tmem30_aminophospholipid_flippase',
        parent = 'transporter',
        resource = 'Almen2009',
        args = {
            'classes': 'TMEM30',
        },
    ),
    af.AnnotDef(
        name = 'tmem16_phospholipid_scramblase',
        parent = 'transporter',
        resource = {'Q6IWH7', 'Q4KMQ2', 'A1A5B4'},
    ),


    # ion channels
    af.AnnotDef(
        name = 'ion_channel',
        resource = '~ion_channel',
        scope = 'generic',
        source = 'composite',
        transmitter = False,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'ion_channel',
        scope = 'generic',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'ion channels',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'ion_channel',
        resource = 'DGIdb',
        scope = 'generic',
        args = {
            'category': 'ION CHANNEL',
        },
        exclude = {
            'Q9ULQ1', 'Q9BRX2', 'Q14289', 'Q9H6F2', 'O14791', 'P14780',
            'F7VJQ1', 'O43768', 'P78509', 'Q8TE54', 'Q96S66', 'Q92508',
            'P78352', 'Q12959', 'Q13976', 'Q9Y4I1', 'P22466', 'P62942',
            'Q9P246', 'P04839', 'Q92915', 'P51790', 'P51793', 'P23327',
            'P0CG08', 'O75534', 'P08133', 'Q9ULM6', 'Q6IQ26', 'P21796',
            'Q86YM7', 'Q9NWR8', 'P68106', 'P54257', 'P28161', 'P56539',
            'Q13387', 'P00367', 'Q9NVV0', 'P80108', 'Q9UMX0', 'P01588',
            'Q9BYP7', 'P62258', 'Q13023', 'P33176', 'Q96PH1', 'Q9P2U7',
            'P34998', 'Q09666', 'Q9UM00', 'Q14393', 'P62879', 'Q9HD26',
            'Q9BV40', 'Q14573', 'P53355', 'Q13127', 'P0DP25', 'Q9HC97',
            'O14958', 'P24387', 'Q96SF2', 'Q9NRX4', 'P78417', 'P58400',
            'Q6XPS3', 'Q6ZUT9', 'Q05513', 'P0DP23', 'Q99959', 'O95833',
            'Q9HDC5', 'P21333', 'Q96PU5', 'Q96MG2', 'P17612', 'P46934',
            'P45880', 'Q96NY7', 'Q8NGH8', 'Q13303', 'Q96BR1', 'Q14643',
            'P23297', 'P84074', 'P51797', 'Q9Y277', 'P25774', 'Q7LC44',
            'Q06787', 'P61328', 'O43448', 'Q9Y696', 'P26678', 'Q14571',
            'P54284', 'Q6ZRF8', 'P06850', 'P01303', 'B7ZAQ6', 'O00141',
            'Q14644', 'Q9BVC6', 'Q13972', 'Q9H4A3', 'O75628', 'Q96D96',
            'Q6TFL4', 'Q8NE86', 'Q9UEU0', 'Q92796', 'Q93034', 'P63165',
            'P02760', 'P49768', 'Q9BSW2', 'P55042', 'P56211', 'Q8NHX9',
            'Q9BYB0', 'Q92913', 'P07550', 'Q8N4C8', 'P05067', 'P35609',
            'Q8NBP7', 'P05771', 'P57796', 'Q92736', 'P29475', 'P19429',
            'Q16623', 'Q8NGH5', 'Q9NZ94', 'Q15413', 'P57727', 'O60733',
            'P30626', 'Q8NGS4', 'P56180', 'P21817', 'Q9Y566', 'Q99653',
            'P42858', 'Q03135', 'P04156', 'O15400', 'Q8N5I3', 'Q9Y6N3',
            'Q8TBE1', 'P62166', 'P58401', 'P01160', 'Q9UBK2', 'Q71RS6',
            'P20936', 'Q9Y6X2', 'P16885', 'Q9P2S2', 'P51798', 'Q8N335',
            'Q9ULB1', 'Q99996', 'P35462', 'Q08499', 'Q8N144', 'Q9BSA9',
            'Q9Y2W7', 'P06756', 'Q8WXH2', 'Q9Y217', 'P30989', 'P0DP24',
            'O75052', 'Q8TEL6', 'P11532', 'Q06413', 'Q16651', 'O14775',
            'Q9HBY8',
        },
    ),
    af.AnnotDef(
        name = 'ion_channel',
        parent = 'ion_channel',
        resource = 'Adhesome',
        scope = 'generic',
        args = {'mainclass': 'Channel'},
    ),  # only 5 channels but is all right
    af.AnnotDef(
        name = 'ion_channel',
        resource = '~ion_channel~HGNC',
        resource_name = 'HGNC',
        scope = 'generic',
    ),
    af.AnnotDef(
        name = 'ion_channel',
        parent = 'ion_channel',
        resource = 'Almen2009',
        scope = 'generic',
        args = {
            'classes': 'Channels',
        },
        exclude = {
            'P51797', 'Q92736', 'Q15413', 'P21817', 'Q14573',
            'Q13520', 'P51793', 'P51798', 'O94778', 'Q8NHX9',
        },
    ),  # all ion channels a few localized only intracellularly
    af.AnnotDef(
        name = 'ion_channel',
        parent = 'ion_channel',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Ion channel',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    # specific subclasses from UniProt
    af.AnnotDef(
        name = 'ligand_gated_channel',
        parent = 'ion_channel',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Ligand-gated ion channel',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'chloride_channel',
        parent = 'ion_channel',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Chloride channel',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'calcium_channel',
        parent = 'ion_channel',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Calcium channel',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'potassium_channel',
        parent = 'ion_channel',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Potassium channel',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'sodium_channel',
        parent = 'ion_channel',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Sodium channel',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    # specific subclasses from Almen 2009
    af.AnnotDef(
        name = 'sperm_associated_voltage_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'catSper-Two-P',
        },
        exclude = {'Q9ULQ1', 'Q8NHX9'},
    ),
    af.AnnotDef(
        name = 'cyclic_nucleotide_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'cycliNucleotideRegulatedChannels',
        },
    ),
    af.AnnotDef(
        name = 'polycystin_calcium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'PKD1',
        },
    ),  # components of calcium channels
    af.AnnotDef(
        name = 'tmem63_osmosensitive_cation',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'TMEM63',
        },
        exclude = {'O94886'},
    ),
    af.AnnotDef(
        name = 'tmem16_calcium_dependent_chloride',
        parent = 'ion_channel',
        resource = {'Q9NQ90'},
    ),
    af.AnnotDef(
        name = 'voltage_gated_calcium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': {'CACNA2D', 'CACNG'},
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_chloride',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'CLC',
        },
    ),
    af.AnnotDef(
        name = 'calcium_activated_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Ca_Activated_Potassium_Channels',
        },
    ),
    af.AnnotDef(
        name = 'atp_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'ATP_gated_ion_channels',
        },
    ),
    af.AnnotDef(
        name = 'calcium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'CalciumChannels',
        },
    ),
    af.AnnotDef(
        name = 'inward_rectifier_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'InwardlyRectifyingKchannel',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'KCNE',
        },
    ),
    af.AnnotDef(
        name = 'calcium_activated_potassium_b',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'KCNMB',
        },
    ),
    af.AnnotDef(
        name = 'ligand_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Ligand_gated_ion_channels',
        },
    ),
    af.AnnotDef(
        name = 'potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Potassium_channels',
        },
    ),
    af.AnnotDef(
        name = 'beta_sodium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'SodiumChannelBeta',
        },
    ),
    af.AnnotDef(
        name = 'sodium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'SodiumChannels',
        },
    ),
    af.AnnotDef(
        name = 'transient_receptor_potential_cation',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'TRP_channels',
        },
    ),  # calcium permeable cation channels activated by receptors
    af.AnnotDef(
        name = 'tweety_chloride',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Tweety',
        },
    ),
    af.AnnotDef(
        name = 'sodium_potassium',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Two-p_K_channels',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated',
        parent = 'ion_channel',
        resource = 'Almen2009',
        args = {
            'classes': 'Voltage_gated_ion_channels',
        },
        exclude = {'Q92736', 'P21817', 'Q14573', 'Q9ULQ1'},
    ),

    # specific subclasses from HGNC
    af.AnnotDef(
        name = 'acid_sensing',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Acid sensing ion channel subunits',
        },
    ),
    af.AnnotDef(
        name = 'calcium_voltage_gated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Calcium voltage-gated channel alpha1 subunits',
                (
                    'Calcium voltage-gated channel auxiliary '
                    'alpha2delta subunits'
                ),
                ' Calcium voltage-gated channel auxiliary beta subunits',
            },
        },
    ), # these look like all being in cell membrane, but maybe we should
       # filter them?
    af.AnnotDef(
        name = 'catsper',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cation channels sperm associated',
        },
    ),
    af.AnnotDef(
        name = 'chloride',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chloride channels, ATP-gated CFTR',
        },
    ),
    af.AnnotDef(
        name = 'cyclic_nucleotide_gated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cyclic nucleotide gated channels',
        },
    ),  # channels gated by various cyclic nucleotides, playing roles in
        # mostly in visual and olfactory signaling
    af.AnnotDef(
        name = 'hydrogen_voltage_gated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Hydrogen voltage gated channels',
        },
    ),
    af.AnnotDef(
        name = 'b_lymphocyte_calcium',
        parent = 'ion_channel',
        resource = {'P11836'},
    ),
    af.AnnotDef(
        name = 'orai_calcium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'ORAI calcium release-activated calcium modulators',
        },
    ),
    af.AnnotDef(
        name = 'calcium_activated_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Potassium calcium-activated channel subfamily '
                'M regulatory beta subunits',
                'Potassium calcium-activated channels',
            },
        },
    ),
    af.AnnotDef(
        name = 'sodium_activated_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Potassium sodium-activated channel subfamily T',
        },
    ),
    af.AnnotDef(
        name = 'two_pore_domain_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Potassium two pore domain channel subfamily K',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_potassium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Potassium voltage-gated channel regulatory subunits',
                'Potassium voltage-gated channel subfamily J',
                'Potassium voltage-gated channels',
            },
        },
        exclude = {
            'O43448', 'Q14722', 'Q9NZI2',
            'Q9Y2W7', 'Q13303', 'Q6PIL6',
        },
    ),
    af.AnnotDef(
        name = 'epithelial_sodium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sodium channels epithelial',
        },
    ),
    af.AnnotDef(
        name = 'sodium_leak',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Sodium leak channels, non selective',
        },
    ),
    af.AnnotDef(
        name = 'voltage_gated_sodium',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': {
                'Sodium voltage-gated channel alpha subunits',
                'Sodium voltage-gated channel beta subunits',
            },
        },
    ),
    af.AnnotDef(
        name = 'transient_receptor_potential_cation',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Transient receptor potential cation channels',
        },
        exclude = {'Q9GZU1'},
    ),
    af.AnnotDef(
        name = 'transmembrane_channel_like',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Transmembrane channel like family',
        },
        exclude = {'Q7Z403', 'Q8IU68'},
    ),
    af.AnnotDef(
        name = 'tweety_chloride',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Tweety family',
        },
    ),  # Ca activated chloride channels
    af.AnnotDef(
        name = 'volume_regulated_anion',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Volume regulated anion channel subunits',
        },
    ),
    af.AnnotDef(
        name = 'zinc_activated',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Zinc activated channels',
        },
    ),


    # ion channel regulator
    af.AnnotDef(
        name = 'ion_channel_regulator',
        resource = '~ion_channel_regulator',
        scope = 'generic',
        source = 'composite',
        transmitter = False,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'chloride',
        parent = 'ion_channel_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Chloride channel accessory',
        },
    ), # in plasma membrane, regulate cholride channels


    # junctions
    # gap junction
    af.AnnotDef(
        name = 'gap_junction',
        scope = 'generic',
        source = 'composite',
        resource = '~gap_junction',
        transmitter = True,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'gap_junction',
        resource = 'GO_Intercell',
        args = {
            'mainclass': 'gap junction',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction',
        resource = 'Ramilowski_location',
        scope = 'generic',
        args = {
            'location': 'gap junction',
        },
        exclude = {'Q69YQ0'},
    ),
    af.AnnotDef(
        name = 'gap_junction',
        resource = 'UniProt_location',
        scope = 'generic',
        args = {
            'location': 'Gap junction',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction',
        resource = 'Almen2009',
        scope = 'generic',
        args = {
            'classes': 'GapJunction',
        },
        exclude = {'A6NN92'},
    ),
    af.AnnotDef(
        name = 'gap_junction',
        scope = 'generic',
        resource = 'HGNC',
        args = {
            'mainclass': 'Gap junction proteins',
        },
    ),
    af.AnnotDef(
        name = 'gap_junction',
        scope = 'generic',
        resource = 'UniProt_keyword',
        args = {
            'keyword': 'Gap junction',
        },
        exclude = {'Q69YQ0', 'P48745'},
    ),
    af.AnnotDef(
        name = 'pannexin',
        parent = 'gap_junction',
        resource = 'HGNC',
        args = {
            'mainclass': 'Pannexins',
        },
    ),  # either half channels or gap junctions


    # tight junction
    af.AnnotDef(
        name = 'tight_junction',
        scope = 'generic',
        source = 'composite',
        resource = '~tight_junction',
        transmitter = True,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'tight_junction',
        resource = 'GO_Intercell',
        scope = 'generic',
        args = {
            'mainclass': 'tight junction',
        },
    ),
    af.AnnotDef(
        name = 'tight_junction',
        resource = 'Ramilowski_location',
        scope = 'generic',
        args = {
            'location': 'tight junction',
        },
        exclude = {
            'Q68EM7', 'Q92974', 'Q9NPG3', 'Q9H8V3', 'O95786', 'Q9NR48',
            'Q92797', 'Q68DX3', 'Q9BUZ4',
        },
    ),
    af.AnnotDef(
        name = 'tight_junction',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Tight junction',
        },
        limit = 'plasma_membrane_transmembrane',
    ),
    af.AnnotDef(
        name = 'tight_junction',
        scope = 'generic',
        resource = 'Zhong2015',
        args = {
            'type': 'tight junction',
        },
    ),
    af.AnnotDef(
        name = 'tight_junction',
        scope = 'generic',
        resource = 'HGNC',
        args = {
            'mainclass': 'Claudins',
        },
    ),
    af.AnnotDef(
        name = 'claudin',
        parent = 'tight_junction',
        resource = 'Almen2009',
        args = {
            'classes': 'Claudin',
        },
    ),


    # adherens junction
    af.AnnotDef(
        name = 'adherens_junction',
        scope = 'generic',
        source = 'composite',
        resource = '~adherens_junction',
        transmitter = False,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'adherens_junction',
        scope = 'generic',
        resource = 'Ramilowski_location',
        args = {
            'location': 'adherens junction',
        },
        exclude = {
            'P35222', 'P35221', 'Q13444', 'Q96TA1', 'Q9UGP4', 'Q9BVG8',
            'A6NIX2', 'Q9HCM4', 'O60331', 'Q9P1Y5', 'Q6IQ23', 'Q8N264',
            'Q8N157',
        },
    ),
    af.AnnotDef(
        name = 'adherens_junction',
        scope = 'generic',
        resource = 'UniProt_location',
        args = {
            'location': 'Adherens junction',
        },
        limit = 'cell_surface',
    ),  # to be checked
    af.AnnotDef(
        name = 'adherens_junction',
        scope = 'generic',
        resource = 'Zhong2015',
        args = {
            'type': 'adherens junction',
        },
    ),


    # desmosome (macula adherens)
    af.AnnotDef(
        name = 'desmosome',
        scope = 'generic',
        source = 'composite',
        resource = '~desmosome',
        transmitter = True,
        receiver = True,
    ),
    af.AnnotDef(
        name = 'desmosomal_cadherin',
        parent = 'desmosome',
        resource = 'HGNC',
        args = {
            'mainclass': 'Desmosomal cadherins',
        },
    ),


    # intracellular protein classes in close relation to intercellular
    # communication
    af.AnnotDef(
        name = 'intracellular_intercellular_related',
        resource = '~intracellular_intercellular_related',
        scope = 'generic',
        source = 'composite',
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'small_molecule_ligand_synthase',
        parent = 'intracellular_intercellular_related',
        scope = 'generic',
        resource = {
            'P09172', # Dopamine beta-hydroxylase (noradrenaline)
            'P35354', # Prostaglandin G/H synthase 2
            'P19113', # Histidine decarboxylase (histamine)
        },
        transmitter = True,
        receiver = False,
    ),
    af.AnnotDef(
        name = 'crumbs_complex',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Crumbs complex',
        },
    ),  # scaffolds and regulators for plasma membrane proteins
    af.AnnotDef(
        name = 'engulfment_motility',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Engulfment and cell motility proteins',
        },
        transmitter = False,
        receiver = True,
    ),  # some intracellular proteins involved in endocytosis
    af.AnnotDef(
        name = 'fbar_actin_dynamics_endocytosis',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'F-BAR domain containing',
        },
    ),  # intracellular proteins, most of them regulate the
        # actin dynamics in endocytosis
    af.AnnotDef(
        name = 'ferm_domain',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'FERM domain containing',
        },
    ),  # intracellular proteins, most of these regulate adhesion and
        # membrane-cytoskeleton interactions; maybe not all related closely
        # to intercellular communication processes
    af.AnnotDef(
        name = 'ferlin',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Ferlin family',
        },
    ),  # intracellular proteins involved in plasma membrane repair
        # and synaptic vesicle fusion
    af.AnnotDef(
        name = 'fermitin',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fermitins',
        },
    ),  # intracellular proteins, peripheral membrane proteins on the
        # cytoplasmic side of the plasma membrane;
        # involved in cell-cell adhesion
    af.AnnotDef(
        name = 'flotillin',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Flotillins',
        },
        transmitter = False,
        receiver = True,
    ),  # intracellular proteins with a role in endocytosis
    af.AnnotDef(
        name = 'arc',
        parent = 'intracellular_intercellular_related',
        resource = {'Q7LC44'},
        transmitter = True,
        receiver = False,
    ),  # intercellular RNA transfer
    af.AnnotDef(
        name = 'interferon_regulator',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Interferon regulatory factors',
        },
    ),  # intracellular proteins mostly transcriptionally 
        # regulating interferons
    af.AnnotDef(
        name = 'cas_scaffold_intracell',
        parent = 'matrix_adhesion',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cas scaffold proteins',
        },
    ),  # intracellular part of cell-matrix (focal) adhesion signaling
    af.AnnotDef(
        name = 'cavin_caveolae',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Cavins',
        },
    ),  # caveolae formation, intercellular
    af.AnnotDef(
        name = 'clathrin_coated_pit',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Clathrin subunits',
        },
    ),  # clathrin coated pit formation, intracellular
    af.AnnotDef(
        name = 'collagen_galactosyltransferase',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Collagen beta(1-O)galactosyltransferases',
        },
    ), # collagen synthesis (in ER)
    af.AnnotDef(
        name = 'junctophilin',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Junctophilins',
        },
    ),  # intracellularily connect the plasma membrane and ER to
        # ensure quick response to membrane potential change
    af.AnnotDef(
        name = 'lims1_adhesion',
        parent = 'intracellular_intercellular_related',
        resource = {'P48059'},
    ),
    af.AnnotDef(
        name = 'maguk_tight_junction',
        parent = 'intracellular_intercellular_related',
        resource = {
            'Q07157', 'Q8N3R9', 'Q9UDY2', 'Q96QZ7', 'Q5T2T1', 'O95049',
        },
    ),  # intracellular scaffolding proteins supporting tight junctions
    af.AnnotDef(
        name = 'parin_adhesion_regulator',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Parvins',
        },
    ),  # intracellular proteins regulating adhesion and integrin signaling
    af.AnnotDef(
        name = 'plakophilin_adhesion_regulator',
        parent = 'intracellular_intercellular_related',
        resource = 'HGNC',
        args = {
            'mainclass': 'Plakophilins',
        },
    ),  # important intracellular parts of cell-cell junctions
    af.AnnotDef(
        name = 'actin_regulation_adhesome',
        parent = 'intracellular_intercellular_related',
        resource = 'Adhesome',
        args = {'mainclass': 'Actin regulation'},
    ),
    af.AnnotDef(
        name = 'adhesion_cytoskeleton_adaptor',
        parent = 'intracellular_intercellular_related',
        resource = 'Adhesome',
        args = {'mainclass': 'Adaptor'},
    ),


    # miscellanous from HGNC -- disabled until decision
    af.AnnotDef(
        name = 'cd_molecule',
        resource = 'HGNC',
        args = {
            'mainclass': 'CD molecules',
        },
        enabled = False,
    ),  # membrane proteins, secreted, receptors, enzymes,
        # adhesion proteins
    af.AnnotDef(
        name = 'c2set_domain',
        resource = 'HGNC',
        args = {
            'mainclass': 'C2-set domain containing',
        },
        enabled = False,
    ),  # these are all plasma membrane proteins, ligands,
        # some receptors and adhesion proteins
    af.AnnotDef(
        name = 'c3_pzp_a2m',
        resource = 'HGNC',
        args = {
            'mainclass': (
                'C3 and PZP like, alpha-2-macroglobulin domain containing'
            ),
        },
        enabled = False,
    ),  # secreted or peripheral on the outer side of plasma membrane
        # enzymes, protease inhibitors, receptors, co-receptors
    af.AnnotDef(
        name = 'tetraspanin_plasma_membrane_regulator',
        resource = 'Almen2009',
        args = {
            'classes': 'Tetraspanin',
        },
        exclude = {'O60635'},
        enabled = False,
    ),  # transmembrane proteins in the plasma membrane, regulate various
        # other proteins such as channels, receptors, adhesion proteins
    af.AnnotDef(
        name = 'bage',
        parent = 'ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'BAGE family',
        },
        enabled = False,
    ),
    af.AnnotDef(
        name = 'cap_ctype_lectin',
        resource = 'HGNC',
        args = {
            'mainclass': 'CAP and C-type lectin domain containing',
        },
        enabled = False,
    ), # these are secreted and affect immune signaling
    af.AnnotDef(
        name = 'cmtm',
        parent = 'receptor_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'CKLF like MARVEL transmembrane domain containing',
        },
        enabled = False,
    ), # transmembrane in plasme membrane; regulate receptor availability
    af.AnnotDef(
        name = 'cacng_ion_channel_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': ' Calcium channel auxiliary gamma subunits',
        },
        enabled = False,
    ),  # transmembrane in plasma membrane; regulate calcium channels and
        # glutamate receptors
    af.AnnotDef(
        name = 'calcium_homeostasis',
        parent = 'ion_channel',
        resource = 'HGNC',
        args = {
            'mainclass': 'Calcium homeostasis modulators',
        },
        enabled = False,
    ),  # taste bud ion and ATP channels
    af.AnnotDef(
        name = 'complement_system_activator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Complement system activation components',
        },
        enabled = False,
    ), # secreted receptors, enzymes and signal transmission proteins
    af.AnnotDef(
        name = 'complement_receptor_and_regulator',
        resource = 'HGNC',
        args = {
            'mainclass': 'Complement system regulators and receptors',
        },
        enabled = False,
    ),  # secreted regulators or membrane bound receptors or inhibitors
        # in the complement system downstream signaling
    af.AnnotDef(
        name = 'fibrinogen_c_domain',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibrinogen C domain containing',
        },
        enabled = False,
    ),  # all are secreted, some of them are ligands, enzymes, other kind of
        # regulators for receptors or adhesion, or ECM proteins
    af.AnnotDef(
        name = 'fibronectin_type_iii',
        resource = 'HGNC',
        args = {
            'mainclass': 'Fibronectin type III domain containing',
        },
        enabled = False,
    ),  # a mixture of plasma membrane transmembrane receptors or adhesion
        # proteins, and also ECM proteins;
        # a few of them are not extracellular at all
        # probably are annotated in other, more specific categories,
        # especially the `Ig-like cell adhesion molecule family`
    af.AnnotDef(
        name = 'immunoglobulin_like',
        resource = 'HGNC',
        args = {
            'mainclass': 'Immunoglobulin like domain containing',
        },
        enabled = False,
    ),  # a mixture of plasma membrane transmembrane receptors or adhesion
        # proteins, and also ECM proteins;
        # a few of them are not extracellular at all
        # probably are annotated in other, more specific categories,
        # especially the `Ig-like cell adhesion molecule family`
    af.AnnotDef(
        name = 'gla_domain',
        resource = 'HGNC',
        args = {
            'mainclass': 'Gla domain containing',
        },
        enabled = False,
    ),  # all secreted, various regulators of blood coagulation, ECM,
        # some enzymes or ligands or regulators of other ligands
    af.AnnotDef(
        name = 'hla',
        parent = 'cell_surface_ligand',
        resource = 'HGNC',
        args = {
            'mainclass': 'Histocompatibility complex',
        },
        enabled = False,
    ),  # histocompatibility antigen complex members for presenting
        # antigens on the cell surface
    af.AnnotDef(
        name = 'vset_domain_containing',
        resource = 'HGNC',
        args = {
            'mainclass': 'V-set domain containing',
        },
        enabled = False,
    ),  # various ligands, receptors and adhesion molecules

)
