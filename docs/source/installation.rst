Installation
============

To install FRMatch via Github: 

.. code-block:: console

   $ git clone https://github.com/BeverlyPeng/frmatch.git
   $ cd frmatch
   
Now we create the environment with all the dependencies required to run FRMatch, including the ability to run a jupyter lab notebook version of the tutorial. This was successfully run on a Macbook pro, Xtools installed, Apple M3 pro chip, Sonoma 14.6 operating system.

Additionally, emacs is installed, to faciliate editing. Feel free to inspect the yaml file (environment.yml). Once the environment is created then we activate the environment. 

.. code-block:: console

   $ conda env create -f environment.yml
   $ conda activate environment

Once the environment is activated.  We need to install the FRMatch package.

.. code-block:: console

   $ pip install .

You should be able to "import frmatch".
