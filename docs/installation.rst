Installing scaffold
==================================

::

    # Assumes you have installed omnium in omnium_env conda env.
    # See omnium installation instructions.
    # https://github.com/markmuetz/omnium/docs/installation.rst
    # Active conda env.
    source activate omnium_env

    git clone https://github.com/markmuetz/cloud_tracking
    git clone https://github.com/markmuetz/scaffold_analysis
    cd cloud_tracking
    pip install -e .
    cd ..
    cd scaffold_analysis
    pip install -e .
    cd ..

Testing installation
====================

::

    # Tell omnium about scaffold
    source activate omnium_env
    export OMNIUM_ANALYSER_PKGS=scaffold
    # If you have the u-an388 data, do this to run realistic analysis tests
    export SCAFFOLD_SUITE_UAN388_DIR=/home/markmuetz/omnium_test_suites/scaffold_test_suite/u-an388

    # Should show current onmium version
    omnium version
    # Should show scaffold analysers
    omnium ls-analysers
    # Will run all tests for scaffold.
    omnium test -a scaffold
