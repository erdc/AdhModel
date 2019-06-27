import os
import param
import logging

from holoviews import Path
from geoviews import Path as GeoPath

from genesis.model import Model
from .mesh import AdhMesh

log = logging.getLogger('AdhModel.adh_model')


class AdhModel(Model):
    """
     Class object to hold all data related to an adh simulation.

     The adhModel object primarily stores, modifies, and outputs data. It contains pass through methods for
     reading, and writing the mesh, hotstart, boundary condition, and result files. Included are methods to
     read individual AdH files or read a suite of model files. There is also a validate method to ensure that
     the data within the model is valid and consistent across internal objects.

     An AdhModel object contains a suite of model parameters and one AdhMesh object. The AdhMesh object
     contains the mesh itself and an AdhSimulation object. The AdHSimulation object stores the boundary
     condition information, the hotstart, and the results - all the data for an individual model run.
     """

    model_name = param.String(
        default="adh_model",
        doc="File name prefix for ADH input files.",
    )
    model_path = param.Foldername(
        default=os.getcwd(),
        doc="Path on disk where ADH files are located.",
    )
    version = param.Number(
        default=5.0,
        bounds=(4.5, None),
        softbounds=(4.5, 5.0),
        doc="Version of AdH executable"
    )
    project_name = param.String(
        default='default',
        doc='Global project name'
    )

    units = param.ObjectSelector(
        default='meters',
        objects=['meters', 'feet']
    )

    path_type = param.ClassSelector(default=GeoPath, class_=Path, is_instance=False, doc="""
             The element type to draw into.""")

    def __init__(self, **params):
        super(AdhModel, self).__init__(**params)
        self.mesh = AdhMesh(crs=self.projection.get_crs())

    @property
    def simulation(self):
        return self.mesh.current_sim

    def read(self, path, project_name='*', crs=None, fmt='nc'):
        """Read in AdH model files as an xarray.Dataset object

        Args:
            path(str, required):
                path to the AdH project files
            project_name(str, optional, default='*'):
                the root name of the AdH project. If not specified then it will be derived
                from the first mesh file found in `path`.
            crs(cartopy.CRS, optional, default=None):
                The projection of the mesh file.
            fmt(stt, optional, default='nc'):
                The format of the file being passed in. Valid options are ['nc', '2dm', '3dm']

        Returns:
            xarray.Dataset variables for the nodes, mesh elements, output datasets and hot-start file datasets
        """
        #TODO look at filename for default format?
        fmts = {
            'nc': self.from_netcdf,
            'ascii': self.from_ascii,
            '2dm': self.from_ascii,
            '3dm': self.from_ascii
        }
        return fmts[fmt](path=path, project_name=project_name, crs=crs)

    def write(self, path, fmt='nc'):
        if fmt != 'nc':
            raise IOError('The only option currently available is nc (netcdf)')
        else:
            # write mesh
            self.write_mesh(file_name=path, fmt=fmt)
            # write hotstart
            self.write_hotstart(file_name=path)
            # write boundary conditions
            self.write_bc(file_name=path, validate=True, fmt='bc')
            # write results
            self.write_results(file_name=path, fmt='nc')

    def from_netcdf(self, *args, **kwargs):
        """Read suite of model files from netcdf file and store data in this model object
        NOTE: Boundary conditions are not stored in netcdf, so they must be read from *.bc file

        Args:
            *args: variable length argument list.
            **kwargs: arbitrary keyword arguments.

        Returns:
            None
        """
        nc_file = os.path.join(f'{kwargs["path"]}', f'{kwargs["project_name"]}.nc')
        # read mesh
        self.read_mesh(*args, **kwargs)
        # read hotstart
        self.read_hotstart(path=nc_file, fmt='nc')
        # read boundary conditions (must be read as ascii) # todo add warning?
        bc_file = os.path.join(f'{kwargs["path"]}', f'{kwargs["project_name"]}.bc')
        self.read_bc(bc_file, fmt='bc')
        # read results
        self.read_results(path=nc_file, fmt='nc')

    def from_ascii(self, *args, **kwargs):
        """Read suite of model files from ASCII file and store data in this model object

        Args:
            *args: variable length argument list.
            **kwargs: arbitrary keyword arguments.

        Returns:
            None
        """
        # set the mesh file name
        mesh_file = os.path.join(f'{kwargs["path"]}', f'{kwargs["project_name"]}.3dm')
        # read mesh
        self.read_mesh(mesh_file, project_name=kwargs['project_name'], crs=kwargs['crs'], fmt='3dm')
        # set hotstart file name
        hot_file = os.path.join(f'{kwargs["path"]}', f'{kwargs["project_name"]}.hot')
        # read hotstart
        self.read_hotstart(path=hot_file, fmt='ascii')
        # read boundary conditions
        bc_file = os.path.join(f'{kwargs["path"]}', f'{kwargs["project_name"]}.bc')
        self.read_bc(bc_file, fmt='bc')
        # read results
        self.read_results(kwargs['path'], project_name=kwargs['project_name'], fmt='ascii')

    def read_mesh(self, *args, **kwargs):
        return self.mesh.read(*args, **kwargs)

    def write_mesh(self, *args, **kwargs):
        return self.mesh.write(*args, **kwargs)

    def read_bc(self, *args, **kwargs):
        return self.simulation.read_bc(*args, **kwargs)

    def write_bc(self, *args, **kwargs):
        return self.simulation.write_bc(*args, **kwargs)

    def read_hotstart(self, *args, **kwargs):
        return self.simulation.read_hotstart(*args, **kwargs)

    def write_hotstart(self, *args, **kwargs):
        return self.simulation.write_hotstart(*args, **kwargs)

    def read_results(self, *args, **kwargs):
        return self.simulation.read_results(*args, **kwargs)

    def write_results(self, *args, **kwargs):
        return self.simulation.write_results(*args, **kwargs)

    def read_result(self, *args, **kwargs):
        return self.simulation.read_result(*args, **kwargs)

    def write_result(self, *args, **kwargs):
        return self.simulation.write_result(*args, **kwargs)

    def validate(self):
        # ensure mesh units and model units match
        if self.units != self.mesh.units:
            log.warning('Model units do not match mesh units')


class AdhModelCoupled(AdhModel):
    """
    Class object for holding AdhModel information that to be used for coupling models.
    """
    def __init__(self, **params):
        super(AdhModelCoupled, self).__init__(**params)

        self.configurationMode = None  # Str: Read previous solution, initialize from previous, read start, new

        # Calculation configuration
        self.numberOfProcessors = None  # Int: Number of processors available for the simulation run

        # Coupling to STWave controls
        self.coupleMode = None  # Str: Couple mode to STWave
        self.coupleValue = None  # Str: Simulation time or number of couples, depending on mode

        # # Output time parameters
        self.timeOutputSeriesMethod = None  # Str: Method for creating the output time series
        self.timeOutputSeriesIncrement = None  # Str: Time increment if manually constructing the timestep series
        self.timeOutputSeriesFilename = None  # Str: Input filename for a custom created series read from file
        self.timeOutputSeries = []  # Constructed output time series

        self.simulationStartTime = None  # Int: Start time of the current simulation
        self.simulationStopTime = None  # Int: Stop time of the current simulation
        self.timeInputSeriesIncrement = None  # Int: Time increment if manually constructing the timestep series
