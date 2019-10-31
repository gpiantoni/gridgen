Grid Generation
===============
Grids are generated from a starting electrode. 
Then subsequent electrodes are added to it.
The order in which the electrodes are added is important, so there are different indexing approaches.

From a Corner
-------------
The most intuitive approach is starting from a corner and then add one value to each column and row.
Note that the numbering is important and need to match the numbering of the actual grid.

North-West
~~~~~~~~~~

.. code-block:: python

    index_corner(2, 4, 'NW')

will give x and y indices following this order

.. csv-table:: 
    :file: ../../tests/data/generators/corner_nw.csv 


North-East
~~~~~~~~~~

.. code-block:: python

    index_corner(2, 4, 'NE')

will give x and y indices following this order

.. csv-table:: 
    :file: ../../tests/data/generators/corner_ne.csv 

South-West
~~~~~~~~~~

.. code-block:: python

    index_corner(2, 4, 'SW')

will give x and y indices following this order

.. csv-table:: 
    :file: ../../tests/data/generators/corner_sw.csv 


South-East
~~~~~~~~~~

.. code-block:: python

    index_corner(2, 4, 'SE')

will give x and y indices following this order

.. csv-table:: 
    :file: ../../tests/data/generators/corner_se.csv 


From the Center
---------------
You can also start from the center.

Up-Down
~~~~~~~
This approach adds one electrode below the reference electrode, then one above the reference electrode, then below again and above again.
Once a column is done, it moves left and right in the same way.

.. code-block:: python

    index_up_down(5, 4)

will give x and y indices following this order

.. csv-table:: 
    :file: ../../tests/data/generators/center_updown.csv 
