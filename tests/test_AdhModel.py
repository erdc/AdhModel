import os
from adhmodel.adh_model import AdhModel

import holoviews.plotting.bokeh
import geoviews.plotting.bokeh

TESTDIR = os.path.dirname(os.path.abspath(__file__))


def test_adhmodel_initialize():
    adhmod = AdhModel(version=5.0)
    assert adhmod.version == 5


def test_adhmodel_read():
    pass


def test_adhmodel_write():
    pass


def test_adhmodel_fromnnetcdf():
    # set test project name
    project_name = 'SanDiego'
    # set project directory
    path = os.path.join(TESTDIR, 'test_files', project_name)
    # instantiate empty adh model object
    adhmod = AdhModel()
    # read model files
    adhmod.from_netcdf(path=path, project_name=project_name, crs=None)
    # todo perform checks on what was loaded


def test_adhmodel_fromascii():
    pass


def test_adhmodel_read_mesh():
    # set test project name
    project_name = 'SanDiego'
    # set project directory
    path = os.path.join(TESTDIR, 'test_files', project_name)
    # instantiate empty adh model object
    adhmod = AdhModel()
    # read the mesh file
    adhmod.read_mesh(path=path, project_name=project_name, crs=None, fmt='nc')
    # todo perform checks on what was loaded


def test_adhmodel_write_mesh():
    pass


def test_adhmodel_read_bc():
    # set test project name
    project_name = 'SanDiego'
    # set project directory
    path = os.path.join(TESTDIR, 'test_files', project_name, f'{project_name}.bc')
    # instantiate empty adh model object
    adhmod = AdhModel()
    # read the mesh file
    adhmod.read_bc(file_name=path, fmt='bc')
    # todo perform checks on what was loaded


def test_adhmodel_write_bc():
    pass


def test_adhmodel_read_hotstart():
    # set test project name
    project_name = 'SanDiego'
    # set project directory
    path = os.path.join(TESTDIR, 'test_files', project_name, f'{project_name}.nc')
    # instantiate empty adh model object
    adhmod = AdhModel()
    # read the mesh file
    adhmod.read_hotstart(path=path, fmt='nc')
    # todo perform checks on what was loaded


def test_adhmodel_write_hotstart():
    pass


def test_adhmodel_read_results():
    pass


def test_adhmodel_write_results():
    pass


def test_adhmodel_read_result():
    pass


def test_adhmodel_write_result():
    pass


def test_adhmodel_validate():
    pass







# def test_load_model():
#     # test_file_path = './SanDiego'
#     # test_file_name = 'SanDiego'
#     # adh_model = AdhModel()
#     # numutils_model, attributes = AdhModel.load_model(path=test_file_path, project_name=test_file_name, netcdf=True)
#     pass  # todo
