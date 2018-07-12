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

Installing lbhem fix
--------------------

There is a problem with the idealised dump files produced by the UM in version 11.0(+):

    I've been trying and failing to load UM dumps in iris, and I think that the cause is that LBHEM (index 17 in the header, from UMDP F03) is being incorrectly set for the idealised UM. At least some fields have LBHEM set to 4, when it should be 3 I think. It is causing iris to not be able to load UM dump files, and I suspect that it only applies to dumps from the idealised UM. In a similar vein, I think LBCODE (index 16) should be set to 9, not 1 as it is. This is using UM vn11.0, and iris vn1.13.0, my suite is u-az437 (based on Adrian Lock's u-aw914). Also, I believe this behaviour has changed since UM vn10.7, as I can load dumps from idealised runs find from that version.

    Is this an issue? I can't see what to do to change this in the UM, and don't have a good idea about the ramifications of such a change anyway. I can patch iris to get round it for the time being.

The fix is as follows:

    Fix lbhem == 4 problem.

    This is code LBHEM (index 17 in UMDP F03), it should equal 3 for
    idealised fields. However, it's been set to 4, which causes iris to
    think that the field is ciruclar in x (wraps around the globe).
    This fix lets you load iris, then set iris.site_configuration key
    'fix_lbhem_eq_4' and this will modify the value of lbhem on load to the
    value fo the key (should be 3). Exposes iris.HAS_LBHEM_FIX to let lib
    know it's there.

To fix, checkout copy of iris 1.13.0 and patch:
DOES NOT WORK ON ARCHER.

::

    source activate omnium_env
    conda uninstall iris
    git clone https://github.com/SciTools/iris
    cd iris
    git checkout v1.13.0
    git checkout -b fix_lbhem_eq_4
    patch -p1 < ../scaffold_analysis/installation/iris-1.13.0_git_diff.fix_lbhem_eq_4.patch
    pip install .

To fix, patch local copy of iris 1.13.0:
(ARCHER HACK).

::

    cd ~/work/anaconda3/envs/omnium_env/lib/python3.6/site-packages/iris/
    patch -p1 < ~/work/scaffold_analysis/installation/iris-1.13.0_git_diff.fix_lbhem_eq_4.lib.patch
    cd ~/work/python3_libs/lib/python3.4/site-packages/Iris-1.13.0-py3.4.egg/iris
    patch -p1 < ~/work/scaffold_analysis/installation/iris-1.13.0_git_diff.fix_lbhem_eq_4.lib.patch
