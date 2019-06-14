import os
import filecmp
import pandas as pd
from pandas.util.testing import assert_frame_equal
from adhmodel.simulation.simulation import BoundaryConditions

import holoviews.plotting.bokeh
import geoviews.plotting.bokeh

TESTDIR = os.path.dirname(os.path.abspath(__file__))


def test_bc_read():
    # bc full file path
    file_name = os.path.join(TESTDIR, 'test_files', 'base_RR_Fine_12052017_1038.bc')
    # instantiate empty boundary conditions object
    bc = BoundaryConditions()
    # read boundary condition file and store into object
    bc.read(file_name)


def test_bc_write(tmpdir):
    # bc full file path
    file_name = os.path.join(TESTDIR, 'test_files', 'base_RR_Fine_12052017_1038.bc')
    # instantiate empty boundary conditions object
    bc = BoundaryConditions()
    # read boundary condition file and store into object
    bc.read(file_name)
    # test directory
    out_file = tmpdir.mkdir('out').join('out.bc')
    # read boundary condition file and store into object
    bc.write(out_file)
    # compare output files
    assert filecmp.cmp(file_name, out_file) is True


def test_bc_add_material():
    pass


def test_bc_get_strings():
    pass


def test_bc_set_defaults_sw2d():
    # instantiate empty boundary condition object
    bc_class = BoundaryConditions()
    # set defaults
    bc_class.set_default_sw2(mesh=None)


def test_bc_validate_sw2d():
    # instantiate empty boundary condition object
    bc_class = BoundaryConditions()
    # set defaults
    bc_class.set_default_sw2(mesh=None)
    # validate parameters
    bc_class.validate_sw2()


def test_bc_check_solution_controls():
    pass


def test_create_time_series():
    # instantiate empty boundary condition object
    bc_class = BoundaryConditions()
    df = pd.DataFrame(data={'X': [0, 120], 'Y': [9.4, 83]})
    bc_class.create_time_series(series_type='SERIES BC', time_series=df, series_id=3)
    assert_frame_equal(bc_class.time_series[3].time_series, df)
