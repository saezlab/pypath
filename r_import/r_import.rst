
.. raw:: html

   <h1>

How to import OmniPath data into R igraph

.. raw:: html

   </h1>

.. code:: r

    source('r_import.r')


.. parsed-literal::

    Loading required package: dplyr
    
    Attaching package: ‘dplyr’
    
    The following objects are masked from ‘package:igraph’:
    
        as_data_frame, groups, union
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    


.. code:: r

    op <- omnipath_graph()

.. code:: r

    op <- omnipath_graph(directed = TRUE)

.. code:: r

    is.directed(op)



.. raw:: html

    TRUE


.. code:: r

    V(op)[1]



.. parsed-literal::

    + 1/6997 vertex, named, from f417565:
    [1] Q13616


.. code:: r

    V(op)[1]$label



.. raw:: html

    'CUL1'


.. code:: r

    E(op)[19]$is_stimulation



.. raw:: html

    1


.. code:: r

    E(op)[19]$sources



.. raw:: html

    <ol>
    	<li><ol class=list-inline>
    	<li>'HPRD'</li>
    	<li>'Laudanna_effects'</li>
    	<li>'DIP'</li>
    	<li>'BioGRID'</li>
    	<li>'SignaLink3'</li>
    	<li>'ACSN'</li>
    	<li>'IntAct'</li>
    	<li>'InnateDB'</li>
    	<li>'CancerCellMap'</li>
    	<li>'Laudanna_sigflow'</li>
    	<li>'Wang'</li>
    </ol>
    </li>
    </ol>



.. code:: r

    E(op)[19]$references



.. raw:: html

    <ol>
    	<li><ol class=list-inline>
    	<li>'18805092'</li>
    	<li>'11961546'</li>
    	<li>'12504026'</li>
    	<li>'11961546'</li>
    	<li>'11956208'</li>
    	<li>'10648623'</li>
    	<li>'12417738'</li>
    	<li>'12609982'</li>
    	<li>'11961546'</li>
    	<li>'20399188'</li>
    	<li>'21765416'</li>
    	<li>'12381738'</li>
    	<li>'11956208'</li>
    	<li>'10648623'</li>
    	<li>'11245432'</li>
    	<li>'16275325'</li>
    	<li>'12684064'</li>
    	<li>'21765416'</li>
    	<li>'12706828'</li>
    	<li>'26725323'</li>
    	<li>'11861641'</li>
    	<li>'16123592'</li>
    	<li>'12215511'</li>
    	<li>'16880526'</li>
    	<li>'24949976'</li>
    	<li>'19250909'</li>
    	<li>'10748083'</li>
    	<li>'23770852'</li>
    	<li>'21169563'</li>
    	<li>'15860010'</li>
    	<li>'26540345'</li>
    	<li>'11384984'</li>
    	<li>'10713156'</li>
    	<li>'22474075'</li>
    	<li>'23535663'</li>
    	<li>'17452440'</li>
    	<li>'15361859'</li>
    	<li>'12417738'</li>
    	<li>'17158585'</li>
    	<li>'22822056'</li>
    	<li>'20638939'</li>
    	<li>'19942853'</li>
    	<li>'12504025'</li>
    	<li>'12481031'</li>
    	<li>'12167173'</li>
    	<li>'22767593'</li>
    	<li>'19933270'</li>
    	<li>'26344197'</li>
    	<li>'12609982'</li>
    	<li>'15659098'</li>
    	<li>'19245792'</li>
    	<li>'12565873'</li>
    	<li>'18805092'</li>
    	<li>'15749712'</li>
    	<li>'15448697'</li>
    	<li>'12840033'</li>
    	<li>'11337588'</li>
    	<li>'11027288'</li>
    	<li>'11483504'</li>
    	<li>'16759355'</li>
    	<li>'11961546'</li>
    	<li>'20832730'</li>
    	<li>'12904573'</li>
    	<li>'10230406'</li>
    	<li>'10230407'</li>
    	<li>'22405651'</li>
    	<li>'18826954'</li>
    	<li>'11359933'</li>
    	<li>'19933270'</li>
    	<li>'10713156'</li>
    	<li>'12706828'</li>
    	<li>'12609982'</li>
    	<li>'11956208'</li>
    	<li>'10648623'</li>
    	<li>'11861641'</li>
    	<li>'10535940'</li>
    	<li>'16759355'</li>
    	<li>'11961546'</li>
    	<li>'10230407'</li>
    	<li>'20832730'</li>
    	<li>'10230406'</li>
    </ol>
    </li>
    </ol>



.. code:: r

    E(op)[19]



.. parsed-literal::

    + 1/42503 edge from f417565 (vertex names):
    [1] P62877->Q13616


.. code:: r

    c(V(op)['P62877']$label, V(op)['Q13616']$label)



.. raw:: html

    <ol class=list-inline>
    	<li>'RBX1'</li>
    	<li>'CUL1'</li>
    </ol>



