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

Testing cylc integration
========================

Tests run with ``u-ax437``, set up to run for 2 days (resubs every day) for S0, S4. Run all tests

