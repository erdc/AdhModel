import os
import platform
import subprocess
import logging
import numpy as np
import xarray as xr
import pandas as pd
import param

from .operation_parameters import OperationParameters
from .iteration_parameters import IterationParameters
from .constituent_properties import ConstituentProperties
from .model_constants import ModelConstants
from .material_properties import MaterialProperties
from .output_control import OutputControl
from .time_control import TimeControl
from .time_series import TimeSeries
from .io import number_of_constituents
from .io import read_bc_file
from .io import write_bc_file
from . import parse_hot_file, parse_dat_file, parse_dat_files, append_array_to_dataset


from genesis.mesh import Simulation

log = logging.getLogger('adhmodel.simulation')


class BoundaryConditions(param.Parameterized):

    boundary_strings = param.DataFrame(
        default=pd.DataFrame(data=[], columns=["CARD", "ID", "ID_0", "ID_1"]),
        doc='Columns for this DataFrame ["CARD", "ID", "ID_0", "ID_1"]. This holds information for the following: '
            'NDS ID ID_0 - ID: id of this node string, ID_0: node id that is part of this string, '
            'EGS ID ID_0 ID_1 - ID: id of this edge string, ID_0, ID_1: first, second node ids of element edge,'
            'MDS ID ID_0 ID_1 - ID: id of this mid string, ID_0, ID_1: first, second element ids,'
            'MTS ID - ID: id of this material string.',
    )
    solution_controls = param.DataFrame(
        default=pd.DataFrame(data=[], columns=["CARD", "CARD_2", "STRING_ID", "XY_ID_0", "XY_ID_1", "XY_ID_2"]),
        doc='Columns for this DataFrame ["CARD", "CARD_2", "STRING_ID", "XY_ID_0", "XY_ID_1", "XY_ID_2"]. '
            'This holds information for the following: '
            'NB DIS S_ID XY_ID - Discharge, S_ID: string id, XY_ID: xy series id for discharge, '
            'NB OVL S_ID XY_ID - Flow, S_ID: string id, XY_ID: xy series for flow, flow per unit area for '
                'material strings and flow per unit length for edge strings, '
            'NB OTW S_ID XY_ID - Water surface elevation, S_ID: string id, XY_ID: xy series for elevation, '
            'NB TRN S_ID C_ID XY_ID - Transport, S_ID: string id, C_ID: constituent id, XY_ID: xy series for '
                'constituent concentration (units dependent of the transport type), '
            'NB SPL S_ID XY_ID - Spillway, S_ID: string id, XY_ID: Series for the percent (%) flow out, '
            'NB OUT O_ID I_ID XY_ID - Flow output inside grid, O_ID: outflow edge id, I_ID: inflow edge id, '
                'XY_ID: series outflow id, '
            'NB TID S_ID - Tidal boundary, S_ID: string id, '
            'DB OVL S_ID X_ID Y_ID - Velocity (2D), S_ID: string id, X_ID: series id for x-velocity, Y_ID: series'
                'id for y-velocity, '
            'DB OVH S_ID X_ID Y_ID DEPTH_ID - Velocity and Depth, S_ID: string id, X_ID: series id for x-velocity, '
                'Y_ID: series id for y-velocity, DEPTH_ID: series id for the depth'
            'DB TRN S_ID C_ID XY_ID - Transport Concentration, S_ID: string id, C_ID: constituent id, XY_ID: '
                'xy series for constituent concentration, '
            'DB LDE S_ID XY_ID - Stationary lid elevation, S_ID: string id, XY_ID: xy series for elevation, '
            'DB LDH S_ID XY_ID - Depth of water under stationary lid, S_ID: string id, XY_ID: xy series for depth, '
            'DB LID S_ID XY_ID - Floating stationary object, S_ID: string id, XY_ID: xy series for the draft of lid, '
            'BR USR S_ID XY_ID - User defined breach displacement, S_ID: string id, XY_ID: xy series for elevation, '
            'OB OF S_ID - Natural Outflow, S_ID: string id, '
            'OFF S_ID - Deactivate string, S_ID: string id, '
            'Not represented here: NB SDR, BR JAI, BR SAS, BR MLM, BR FRO, BR BRC, BR VTG, BR FER, WER, WRS, FLP, '
            'FGT, SLUICE, SLS.',
    )
    stage_discharge_boundary = param.DataFrame(
        default=pd.DataFrame(data=[], columns=['CARD', 'S_ID', 'COEF_A', 'COEF_B', 'COEF_C', 'COEF_D', 'COEF_E']),
        doc='Columns [CARD, S_ID, COEF_A, COEF_B, COEF_C, COEF_D, COEF_E], '
            'NB SDR: Stage discharge boundary. S_ID: String id, COEF_A: coeficient A, COEF_B: coeficient B, '
            'COEF_C: coeficient C, COEF_D: coeficient D, COEF_E: coeficient E.',
    )
    friction_controls = param.DataFrame(
        default=pd.DataFrame(data=[], columns=["CARD", "CARD_2", "STRING_ID", "REAL_01", "REAL_02",
            "REAL_03", "REAL_04", "REAL_05"]),
        doc='Columns for this DataFrame ["CARD", "CARD_2", "STRING_ID", "REAL_01", "REAL_02", '
            '"REAL_03", "REAL_04", REAL_05"]. '
            'This holds information for the following: '
            "FR MNG S_ID MAN_N - Manning's N, S_ID: string id, MAN_N: manning's n, "
            "FR MNC S_ID MAN_N - Manning's Equation, S_ID: string id, MAN_N: manning's n, "
            "FR ERH S_ID ROUGH_HEIGHT - Equivalent roughness height, S_ID: string id, ROUGH_HEIGHT: roughness height, "
            "FR SAV S_ID ROUGH_HEIGHT STEM_HEIGHT - Submerged aquatic vegetation, S_ID: string id, ROUGH_HEIGHT: "
                "roughness height of the SAV canopy (k), STEM_HEIGHT: undeflected stem height, "
            "FR URV S_ID ROUGH_HEIGHT STEM_DIAMETER STEM_DENS - Un-submerged rigid vegetation, S_ID: string id, "
                "ROUGH_HEIGHT: bed roughness height, STEM_DIAMETER: avg stem diameter, STEM_DENS: avg stem density, "
            "FR EDO S_ID ROUGH_HEIGHT OBSTRUCT_DIAMETER OBSTRUCT_HEIGHT - Equivalent drag obstructions, "
                "S_ID: string id, ROUGH_HEIGHT: bed roughness height, OBSTRUCT_DIAMETER: obstruction diameter, "
                "OBSTRUCT_HEIGHT: obstruction height, "
            "FR ICE S_ID ICE_THICKNESS ICE_DENSITY ICE_MOVEMENT - Ice thickness, S_ID: string id, "
                "ICE_THICKNESS: ice thickness, ICE_DENSITY: ice density, ICE_MOVEMENT: ice movement, "
            "FR IRH S_ID ICE_ROUGHNESS - Ice roughness, S_ID: string id, ICE_ROUGHNESS: ice roughness height, "
            "FR BRH S_ID BED_ROUGHNESS - Ice bed roughness height, S_ID: string id, BED_ROUGHNESS: bed "
                "roughness height, "
            "FR SDK S_ID HEIGHT - Submerged dike, S_ID: string id, HEIGHT: height of the dike, "
            "FR BRD S_ID ELEV THICKNESS - Bridge deck, S_ID: string id, ELEV: elevation of bridge deck, "
                "thickness of bridge deck",
    )
    breach_controls = param.DataFrame(
        default=pd.DataFrame(data=[], columns=['CARD', 'CARD_1', 'C_00', 'C_01', 'C_02', 'C_03', 'C_04',
                                               'C_05', 'C_06', 'C_07', 'C_08', 'C_09']),
        doc='[CARD][CARD_1][C_00][C_01]...[C_09]'
            'BR JAI S_ID BREACH_SEC WIDTH MIN_ELEV CREST_ELEV FAIL_TIME FURTHEST CLOSEST - Breach johnson and illes, '
                'S_ID: string id, BREACH_SEC: breach section, WIDTH: width of main breach, MIN_ELEV: min breach elev, '
                'CREST_ELEV: Dam/levee crest elevation, FAIL_TIME: breach failure time, FURTHEST: side slope node '
                'furthest from breach, CLOSEST: side slope node closest to breach, '
            'BR SAS S_ID BREACH_SEC WIDTH MIN_ELEV CREST_ELEV FAIL_TIME FURTHEST CLOSEST - Breach singh and snorrason, '
                'S_ID: string id, BREACH_SEC: breach section, WIDTH: width of main breach, MIN_ELEV: min breach elev, '
                'CREST_ELEV: Dam/levee crest elevation, FAIL_TIME: breach failure time, FURTHEST: side slope node '
                'furthest from breach, CLOSEST: side slope node closest to breach, '
            'BR MLM MAX_WATER_DEPTH RESERVOIR_VOL MIN_ELEV CREST_ELEV - Breach macdonald and landgridge-monopolis, '
                'MAX_DEPTH: Max water depth above breach bottom, RESERVOIR_VOL: reservoir volume, MIN_ELEV: minimum '
                'breach elevation, CREST_ELEV: Dam/levee crest elevation, '
            'BR FRO S_ID BREACH_SEC WIDTH MIN_ELEV CREST_ELEV FAIL_TIME EXPO FURTHEST CLOSEST - Breach froelich, '
                'S_ID: string id, BREACH_SEC: breach section, WIDTH: width of main breach, MIN_ELEV: min '
                'breach elev, CREST_ELEV: Dam/levee crest elevation, FAIL_TIME: breach failure time, EXPO: exponent '
                'for main breach, FURTHEST: side slope node furthest from breach, CLOSEST: side slope node closest '
                'to breach, '
            'BR BRC S_ID BREACH_SEC WIDTH MIN_ELEV CREST_ELEV FAIL_TIME FURTHEST CLOSEST - Breach bureau of '
                'reclamation, S_ID: string id, BREACH_SEC: breach section, WIDTH: width of main breach, MIN_ELEV: min '
                'breach elev, CREST_ELEV: Dam/levee crest elevation, FAIL_TIME: breach failure time, FURTHEST: side '
                'slope node furthest from breach, CLOSEST: side slope node closest to breach, '
            'BR VTG S_ID BREACH_SEC WIDTH MIN_ELEV CREST_ELEV EROSION_CHAR FAIL_TIME FURTHEST CLOSEST - Breach '
                'von thun and gillette, S_ID: string id, BREACH_SEC: breach section, WIDTH: width of main breach, '
                'MIN_ELEV: min breach elev, CREST_ELEV: Dam/levee crest elevation, EROSION_CHAR: erosion character, '
                'FAIL_TIME: breach failure time, FURTHEST: side slope node furthest from breach, CLOSEST: side '
                'slope node closest to breach, '
            'BR FER S_ID BREACH_SEC WIDTH MIN_ELEV CREST_ELEV ENG_CHAR FAIL_TIME FURTHEST CLOSEST - Breach '
                'fed energy regulatory disp, S_ID: string id, BREACH_SEC: breach section, WIDTH: width of main breach, '
                'MIN_ELEV: min breach elev, CREST_ELEV: Dam/levee crest elevation, ENG_CHAR: engineering character, '
                'FAIL_TIME: breach failure time, FURTHEST: side slope node furthest from breach, CLOSEST: side '
                'slope node closest to breach, ',
    )
    weirs = param.DataFrame(
        default=pd.DataFrame(data=[], columns=['WRS NUMBER', 'S_UPSTREAM', 'S_DOWNSTREAM', 'WS_UPSTREAM',
                                               'WS_DOWNSTREAM', 'LENGTH', 'CREST_ELEV', 'HEIGHT']),
        doc='WRS NUMBER S_UPSTREAM S_DOWNSTREAM WS_UPSTREAM WS_DOWNSTREAM LENGTH CREST_ELEV HEIGHT - Weir Parameters, '
                'NUMBER: Weir number, S_UPSTREAM: string upstream of weir, S_DOWNSTREAM: string downstream of weir, '
                'NUMBER: weir number, WS_UPSTREAM: weir string on upstream, WS_DOWNSTREAM: weir string on downstream, '
                'LENGTH: length of weir, CREST_ELEV: crest elevation, HEIGHT: weir height',
    )
    flap_gates = param.DataFrame(
        default=pd.DataFrame(data=[], columns=['FGT NUMBER', 'USER', 'S_UPSTREAM', 'S_DOWNSTREAM', 'FS_UPSTREAM',
                                               'FS_DOWNSTREAM', 'COEF_A', 'COEF_B', 'COEF_C', 'COEF_D', 'COEF_E']),
        doc='FGT NUMBER USER S_UPSTREAM S_DOWNSTREAM FS_UPSTREAM FS_DOWNSTREAM COEF_A COEF_B COEF_C COEF_D COEF_E '
                'COEF_F LENGTH - Flap gate parameters, NUMBER: flap gate number, USER: user defined parameters, '
                'S_UPSTREAM: string upstream of flap, S_DOWNSTREAM: string downstream of flap, FS_UPSTREAM: flap '
                'string on the upstream, FS_DOWNSTREAM: flap string on the downstream, COEF_A: coeficient A, '
                'COEF_B: coeficient B, COEF_C: coeficient C, COEF_D: coeficient D, COEF_E: coeficient E, '
                'COEF_F: coeficient F, LENGTH: length of flap gate',
    )
    sluice_gates = param.DataFrame(
        default=pd.DataFrame(data=[], columns=['SLS NUMBER', 'S_UPSTREAM', 'S_DOWNSTREAM', 'SS_UPSTREAM',
                                               'SS_DOWNSTREAM', 'LENGTH', 'TS_OPENING']),
        doc='SLS NUMBER S_UPSTREAM S_DOWNSTREAM SS_UPSTREAM SS_DOWNSTREAM LENGTH TS_OPENING - Sluice gate parameters, '
                'NUMBER: sluice gate number, S_UPSTREAM: string upstream of sluice gate, S_DOWNSTREAM: string '
                'downstream of sluice gate, SS_UPSTREAM: sluice string on upstream, SS_DOWNSTREAM: sluice string on '
                'downstream, LENGTH: length of sluice gate, TS_OPENING: time series defining the sluice gate opening',
    )

    def __init__(self, **params):
        """
        Initializes and sets default values for attributes in the adhModel object by calling subfunctions for
        initialization based on the version number. Default version is 4.5.


        Initilize an AdH version 5.0 object

        Initializes and sets default values for attributes in the adhModel object, including path controls,
        operating parameters, iteration parameters, mesh refinement parameters, friction parameters,
        solution parameters, time parameters, input and output controls, input and output time parameters,
        mesh configuration, depth data, and linked data inside the adhModel object.

        NOTE: Adh data types that are not yet included:
            WNDLIB: List of wind library parameters
            FLI: Maximum number of linear iterations. Advance after obtaining
            FNI: Number of nonlinear iterations per timestep.  Advance after obtaining.
            RTL: Runge-kutta tolerance for reactive constituents
            SST: Quasi-unsteady tolerance
            BEDLAYERS: Total number of bed layers
            RAD: Dirichlet - short wave radiation and dew point
            SRC: NB SOURCE - rain or evap from an area
            INS: Ice strings
            SXX: Radiation stress tensor, Sxx
            SXY: Radiation stress tensor, Sxy
            SYY: Radiation stress tensor, Syy

        Parameters
        ----------
        self : obj
            The adhModel object for this simulation.

        Returns
        -------
        None
            Sets and modifies boundary condition attributes
        """
        super(BoundaryConditions, self).__init__(**params)

        self.time_control = TimeControl()
        self.output_control = OutputControl()
        self.operation_parameters = OperationParameters()
        self.iteration_parameters = IterationParameters()
        self.constituent_properties = ConstituentProperties()
        self.model_constants = ModelConstants()
        self.material_properties = {1: MaterialProperties()}
        self.time_series = {}

    def read(self, file_name, fmt='bc'):
        """
        Function to read in the boundary condition file into the adhModel object.

        Reads an AdH formatted boundary condition file (*.bc) and stores the parameter data for each variable in the
        appropriate space of the simulation object (self) .

        Parameters
        ----------
        file_name : str
            the full file path for the file that is to read in.
        fmt : str
            format of the file. Options: 'bc'

        Returns
        -------
        None
            Saves the results to the appropriate object in adhModel.

        """
        if fmt == 'bc':
            read_bc_file(file_name, self)
        else:
            raise IOError('Boundary condition file format not recognized. Current options are: "bc"')

    def write(self, file_name, validate=True, fmt='bc'):
        """
        Writes the boundary condition file for the o_adhModel object.

        Writes an AdH formatted boundary condition file (*.bc) using data stored in the o_adhModel object file (does
        not include constituents, sediment, winds, or waves).

        Parameters
        ----------
        file_name : str
            the full file path written for the boundary condition data.
        fmt : str
            format of the file. Options: 'bc'
        validate : bool
            Extra validation variable to make sure the bc file is executed by AdH.

        Returns
        -------
        None
            Writes out boundary condition data to specified file and path.

        """
        if fmt == 'bc':
            write_bc_file(file_name, self, validate=validate)
        else:
            raise IOError('Boundary condition file format not recognized. Current options are: "bc"')

    def add_material(self, material, number=None):
        """
        Set a new material or replace an existing one

        :param material: MaterialProperty()
        :param number: Integer
        :return:

        """
        if number is None:
            number = max(self.material_properties, key=int) + 1

        self.material_properties[number] = material

    def get_strings(self, string_type='EGS'):
        """

        :param string_type: String card. Options are 'EGS', 'NDS', 'MTS', 'MDS'
        :return: dataframe of this string type
        """

        df = self.boundary_strings.loc[self.boundary_strings['CARD'] == string_type]

        return df

    def set_default_sw2(self, mesh=None):
        """
        Method to set default boundary conditions for a d2 AdH simulation. Contains everything required to run a
        simplified simulation. Takes an empty AdH simulation object and fills it with default values.

        Note: only applicable for v5.0

        Parameters
        ----------
        self : AdhModel
            Am empty AdH simulation object.

        Returns
        -------
        None
            Sets operating conditions in the adhModel with default values.


        """
        # For reference:
        # String 1 - material string
        # String 2 - edge string, inflow
        # String 3 - edge string, outflow

        # Solution controls
        self.solution_controls = self.solution_controls.append(
            {'CARD': 'NB', 'CARD_2': 'OVL', 'STRING_ID': 2, 'XY_ID_0': 1, 'XY_ID_1': np.NaN, 'XY_ID_2': np.NaN},
            ignore_index=True)

        # Natural - tailwater elevation
        self.solution_controls = self.solution_controls.append(
            {'CARD': 'NB', 'CARD_2': 'OTW', 'STRING_ID': 2, 'XY_ID_0': 1, 'XY_ID_1': np.NaN, 'XY_ID_2': np.NaN},
            ignore_index=True)

        # Time parameters
        self.time_control.start_time = 0.0  # Start time of the simulation
        self.time_control.start_time_units = '0 - seconds'
        self.time_control.end_time = 3600  # End time of the simulation
        self.time_control.end_time_units = '0 - seconds'

        # Boundary string parameters
        # self.mts.append([1, 1])           # Material strings
        self.boundary_strings = self.boundary_strings.append(
            {'CARD': 'MTS', 'ID': 1, 'ID_0': 1, 'ID_1': np.NaN},
            ignore_index=True)

        # if a mesh was provided
        if mesh is not None:
            # if the mesh already has boundary nodes determined
            if mesh.boundaryNodes is not None:
                # todo add this back in after integration with genesis mesh
                # grab the first 3 nodes and create a string

                # grab the middle 3 nodes and create a string
                raise RuntimeError('utilizing a mesh for setting defaults is not currently available. ')

            else:
                print(
                    'Mesh does not contain boundary node information, setting arbitrary nodes into boundary condition strings')
                self.boundary_strings = self.boundary_strings.append(
                    {'CARD': 'EGS', 'ID': 1, 'ID_0': 2, 'ID_1': 2},
                    ignore_index=True)
                self.boundary_strings = self.boundary_strings.append(
                    {'CARD': 'EGS', 'ID': 4, 'ID_0': 5, 'ID_1': 3},
                    ignore_index=True)
        else:
            print('No mesh was provided, setting arbitrary nodes into boundary condition strings')

            self.boundary_strings = self.boundary_strings.append(
                {'CARD': 'EGS', 'ID': 1, 'ID_0': 2, 'ID_1': 2},
                ignore_index=True)
            self.boundary_strings = self.boundary_strings.append(
                {'CARD': 'EGS', 'ID': 4, 'ID_0': 5, 'ID_1': 3},
                ignore_index=True)

        # Time series object for discharge
        temp_data = pd.DataFrame(data=[[0.0, 100.0], [999999, 100.0]], columns=['X', 'Y'])
        temp_series = TimeSeries(series_id=1, series_type='SERIES BC', time_series=temp_data)
        self.time_series[1] = temp_series

        # Time series object for tailwater
        temp_data = pd.DataFrame(data=[[0.0, 20], [999999, 20]], columns=['X', 'Y'])
        temp_series = TimeSeries(series_id=2, series_type='SERIES BC', time_series=temp_data)
        self.time_series[2] = temp_series

        # Time series object for time step
        temp_data = pd.DataFrame(data=[[0.0, 300], [999999, 300]], columns=['X', 'Y'])
        temp_series = TimeSeries(series_id=3, series_type='SERIES DT', time_series=temp_data)
        self.time_series[3] = temp_series

        # Time series object for output series
        temp_data = pd.DataFrame(data=[[0.0, 999999, 300, 0]],
                                 columns=['START_TIME', 'END_TIME', 'TIME_STEP_SIZE', 'UNITS'])
        temp_series = TimeSeries(series_id=4, series_type='SERIES AWRITE', time_series=temp_data)
        self.time_series[4] = temp_series

    def validate_sw2(self):
        """
        Method to validate AdH Parameters to ensure that the required parameters are included and are reasonable.

        Parameters
        ----------

        Returns
        -------
        None
            Has the potential to stop the program with multiple exceptions. Prints a success statement upon completion.
        """
        # ### Operating parameters ###
        # check preconditioner
        if self.operation_parameters.preconditioner_type != 1:
            log.warning('Preconditioner is not set to the recommended default value of 1. Proceed with caution.')

        # check blocks
        if self.operation_parameters.blocks_per_processor != 1:
            log.warning('Number of blocks is not set to the recommended initial value of 1. Proceed with caution.')

        # ### Iteration Parameters ###
        # check tolerance is reasonable
        if self.iteration_parameters.non_linear_residual_tolerance > 0.5:
            log.warning('Iteration tolerance is unrealistic.')

        # check tolerance is reasonable
        if self.iteration_parameters.non_linear_incremental_tolerance > 0.5:
            log.warning('Iteration tolerance is unrealistic.')

        # ### Material parameters ###
        # check density
        if not (999.0 < self.model_constants.density < 1001.0) and not (
                1.93 < self.model_constants.density < 1.95):
            log.warning('Density should be specified in SI (1000 kg/cm) or English (1.94 slugs/cft). '
                       'Syntax: MP RHO 1000 ')

        # check manning's n constant
        if not self.model_constants.mannings_unit_constant == 1.0 and not (
                1.4 < self.model_constants.mannings_unit_constant < 1.5):
            log.warning('Mannings unit constant should be specified in SI (1.0) or English (1.486). '
                       'Syntax: MP MUC 1.0 ')

        # check material strings
        max_material = max(self.material_properties, key=int)
        if max_material != len(self.material_properties):
            raise IOError('*** Material strings not set properly. Material numbering must be sequential ***')

        # check estimated eddy viscosity
        for key, material in self.material_properties.items():
            reasonable_bounds = material.param.estimated_eddy_viscosity_weighting_factor.softbounds
            if not reasonable_bounds[0] <= material.estimated_eddy_viscosity_weighting_factor <= reasonable_bounds[1]:
                raise log.warning('*** Estimated eddy viscosity is outside of suggested range.***')

            # check transport properties
            if self.operation_parameters.transport == 0:
                if material.transport_properties:
                    raise log.warning(
                        '*** Transport properties given for Material {}, but OP TRN is set to 0'.format(key))
            elif len(material.transport_properties) != self.operation_parameters.transport:
                raise log.warning('*** Missing transport properties for some constituents in Material {}'.format(key))

        # ### Time series ###
        # ensure designated time step series is correct
        if self.time_control.time_step_option == 'Time step series (SERIES DT)':
            if not self.time_series[
                       self.time_control.max_time_step_size_time_series].series_type == 'SERIES DT':
                raise IOError('*** SERIES DT number does not match specification in time control ***')

        # get a list of the series types in the simulation
        type_list = [value.series_type for key, value in self.time_series.items()]
        # check time series objects
        if 'SERIES BC' not in type_list:
            raise IOError('*** No SERIES BC found in the simulation ***')
        if 'SERIES AWRITE' not in type_list and 'SERIES WRITE' not in type_list and 'SERIES OS' not in type_list:
            raise IOError('*** No output series found in the simulation ***')

        # check_transport(self):
        num_trn = number_of_constituents(self)
        op_trn = self.operation_parameters.transport
        if num_trn != op_trn:
            raise IOError('Invalid number of constituents specified. OP TRN {} does not match '
                               'the {} constituents found in the simulation.'.format(op_trn, num_trn))

        errors = self.check_solution_controls()

        # todo: check for boundary condition statements (string/series linkage)
        # this list is HUGE @@ update this

        result = True  # todo collect warnings and return
        return result

    def check_solution_controls(self):
        """
        Checks specified boundary conditions to make sure the material, node, edge, mid string exists and that any
        assigned time series exists
        """
        errors = []
        # create dict of boundary string ids
        df2 = self.solution_controls
        if df2 is None:
            return
        df = self.boundary_strings
        if df is None:
            errors.append('No boundary strings defined.')
            return errors
        if 'ID' not in df:
            errors.append('Invalid definition of boundary strings. No "ID" column defined.')
            return errors
        if 'STRING_ID' not in df2:
            errors.append('Invalid definition of solution controls. No "STRING_ID" column defined.')
            return errors

        bs_set = set()
        for index, row in df.iterrows():
            bs_set.add(row['ID'])
        sid_set = set()
        for index, row in df2.iterrows():
            sid_set.add(row['STRING_ID'])
        for sid in sid_set:
            if sid not in bs_set:
                errors.append('Invalid solution control specified with "STRING ID" {}. This string id is not '
                                   'defined.'.format(sid))
        return errors

    def create_time_series(self, series_type, time_series, units=0, series_id=None, **kwargs):
        """Method to make AdH time series object from an array of time and data.
        Parameters
        ----------
        series_type: str
            Type of time series. Supported types in ADH 5.0 are 'SERIES AWRITE', 'SERIES WRITE', 'SERIES BC',
            'SERIES DT', 'SERIES WIND', 'SERIES WAVE'
        time_series: dataframe
            Time series. Column labels ['X','Y']
        units: int
            time units. Options include:
                0 - seconds
                1 - minutes
                2 - hours
                3 - days
                4 - weeks
        series_id: int
            ID of this xy series
        kwargs
        ------
        output_units: int
            output time units. Options include:
                0 - seconds
                1 - minutes
                2 - hours
                3 - days
                4 - weeks
        x_location: float
            X location for SERIES WIND and SERIES WAVE
        y_location: float
            Y location for SERIES WIND and SERIES WAVE
        Returns
        -------
        None
           Sets self.time_series{key: value}, a time series object inside of adh_model containing all time series information.
        """
        # if no series id was given, append to the end of the existing list of time series
        if series_id is None:
            series_id = len(self.time_series) + 1
        # check to see if series id already exist
        elif series_id in self.time_series.keys():
            log.warning('Overwriting existing time series with same series_id')
        # create new series
        ts = TimeSeries(series_id=series_id, series_type=series_type, time_series=time_series,
                        units=units, **kwargs)
        # set series into model object
        self.time_series[series_id] = ts


class AdhSimulation(Simulation):
    boundary_conditions = param.ClassSelector(default=BoundaryConditions(), class_=BoundaryConditions)
    results = param.ClassSelector(default=xr.Dataset(), class_=xr.Dataset)
    hotstart = param.ClassSelector(default=xr.Dataset(), class_=xr.Dataset)

    def __init__(self, **params):
        super(Simulation, self).__init__(**params)

    def read_bc(self, *args, **kwargs):
        self.boundary_conditions.read(*args, **kwargs)

    def write_bc(self, *args, **kwargs):
        self.boundary_conditions.write(*args, **kwargs)

    def read_hotstart(self, path, fmt='nc'):
        if fmt == 'nc':
            model_xr = xr.open_dataset(path)
            self.hotstart = self._subset_xarray_vars_by_coords(model_xr, ('init_time', 'nodes_ids'))
        elif fmt == 'ascii':
            try:
                self.hotstart = parse_hot_file(path)
            except ValueError as e:
                log.warning(e)

    def write_hotstart(self, file_name):
        """
        Writes the ADH hot start file from xarray

        Args:
            file_name: file name for *.hot file
        """
        with open(file_name, 'w') as hotstart_file:
            for varname, ht in self.hotstart.data_vars.items():
                header = ht.attrs.copy()
                header['OBJTYPE'] = f'"{header["OBJTYPE"]}"'
                header.pop('DIM', None)
                for k, v in header.items():
                    hotstart_file.write(f'{k} {v}\n')
                hotstart_file.write(f'NAME "{ht.name}"\n')
                hotstart_file.write('TS 0 0\n')
                hotstart_file.write(ht.squeeze().to_pandas().to_csv(sep=' ', index=False, header=False,
                                                                    line_terminator='\n'))
                hotstart_file.write('ENDDS\n')

    def read_results(self, path, project_name='*', fmt='nc', **kwargs):
        if fmt == 'nc':
            model_xr = xr.open_dataset(path)
            self.results = self._subset_xarray_vars_by_coords(model_xr, ('times', 'nodes_ids'))
        elif fmt == 'ascii':
            dat_path = os.path.join(path, project_name + '_*.dat')
            if 'file_list' in kwargs:
                file_list = kwargs['file_list']
            else:
                file_list = None
            self.results = parse_dat_files(dat_path, project_name=project_name, file_list=file_list)

    def _subset_xarray_vars_by_coords(self, dataset, coords):
        # determine what attributes are results
        coords_set = set(coords)
        subset_vars = []
        for var in dataset.data_vars:
            dims_set = set(dataset[var].dims)
            if coords_set.issubset(dims_set):
                subset_vars.append(var)
        return dataset[subset_vars]

    def write_results(self, file_name, fmt='nc'):
        if fmt == 'nc':
            self.results.to_netcdf(path=file_name)

    def read_result(self, file_name):
            array = parse_dat_file(file_name)
            append_array_to_dataset(self.results, array)

    def write_result(self, file_name, result_var_name):
        pass

    def solve_adh_model(self, path, filename, hpc_compute=False, architecture='auto'):
        """
        Method to execute AdH based on the information stored in o_adhModel. This method includes capability for
        solving using 32 or 64 bit machines and/or multiprocessor.

        Checks whether simulation is being solved on local PC or HPC. As of this writing of 5/23/17, HPC is not supported,
        and will raise an error stating a message to similar effect. Sets the directory to path, then checks the 
        processors. After checking that operations are within parameters, runs the executables containing the simulation
        solves. Outputs results to the specified path and file. Reports whether the solve ran successfully, and 
        catches most run failures. However, some can get by the system's checks. 

        Parameters
        ----------
        path : str
            Parameter for determining the end location for the file of the solved AdHModel data. 
        filename : str
            Parameter for naming the end file containing the solved adhModel data. 
        hpc_compute : bool
            Keeps track of whether the simulation is being solved on multiprocessor (HPC) or local PC compute.
            If set to True, it is a HPC compute;  if set to False, it is a local PC compute.
        architecture : str
                Used to determine architecture
                auto: Detect platform architecture and run the correct module
                32bit: Run the 32bit executable
                64bit: Run the 64bit executable

        Returns
        -------
        None.
            Writes standard AdH output files (dependent upon simulation parameters).
            'adh.out' is the file containing AdH screen output information.

        Notes
        -----
        Executable file names are hard coded (same as executables on the website) and are expected to
        be in the /adh/ folder [1]_

        Will catch most AdH failed runs, but not all (e.g. won't catch failed run due to missing files) [2]_

        Temporarily will not run HPC uses (05/23/17) [3]_

        Cannot run multi-processor solves on local PC with 32 bit architecture: multi-processing has to be done on
        a local PC with 64 bit architecture as of (05/23/17) [4]_.

        """
        if not (architecture == 'auto' or architecture == '32bit' or architecture == '64bit'):
            raise IOError('AdH solution architecture is not understood.')

        # ### Get information about the current machine ###
        s_operatingSystem = platform.system()

        if architecture == 'auto':
            architecture = platform.architecture()[0]

        ### Get the path to the executables ###
        # Get the path the the current subdirectory
        s_repositorySubpackagePath = os.path.dirname(input_adh_mesh.__file__)

        # Move up one directory by splitting and reconstituting the path
        sl_repositorySubpackagePathSplit = s_repositorySubpackagePath.split('\\')

        s_repositoryPath = sl_repositorySubpackagePathSplit[0] + '\\'
        for i_entry in range(1, len(sl_repositorySubpackagePathSplit) - 1, 1):
            s_repositoryPath = os.path.join(s_repositoryPath, sl_repositorySubpackagePathSplit[i_entry])

        ### Solve the simulation ###
        # If solving on local machine
        if hpc_compute == False:
            # Store plans current directory
            s_curdir = os.getcwd()

            # Change the current working directory
            os.chdir(path)

            # Call the solve routine
            if self.numberOfProcessors is not None and self.numberOfProcessors > 1:
                # Check version
                if self.version == 5.0:
                    raise IOError('Multiprocessor not yet available for AdH v5')

                # Multiprocessing in enabled.  Strip the output across all available cores
                # Set executable paths based on machine architecture
                if s_operatingSystem == 'Windows' and architecture == '64bit':
                    # Construct the path to the windows multicore executables
                    s_executableFullPathPre = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                           'pre_adh_V4.5-MUL.exe')
                    s_executableFullPathAdh = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                           'adh_V4.5-MUL.exe')

                elif s_operatingSystem != 'Windows':
                    # Multi-processing is not supported on non-windows platforms
                    raise NotImplementedError('AdH multiprocessor not available for non-Windows systems')

                else:
                    # Multi-processing is not supporting on 32bit systems
                    raise IOError('AdH multiprocessor not available for 32bit')

                # Execute pre_adh
                # Open the log file
                o_preAdHLog = open('pre_adh.out', 'w')

                # Call the process
                o_preAdHProcess = subprocess.Popen(['mpiexec', '-n', '1', s_executableFullPathPre, filename],
                                                   shell=False, stderr=subprocess.PIPE, stdout=o_preAdHLog)

                # Flush and close the log file
                o_preAdHLog.flush()
                o_preAdHLog.close()

                # Get the pre-AdH errors
                o_error = o_preAdHProcess.communicate()

                # Request the status code from the process
                i_status = o_preAdHProcess.returncode

                # Throw an error if the simulation did not complete successfully
                if i_status != 0:
                    # Error status received.
                    raise IOError('PreAdH failed to validate ' + filename + ' simulation.')

                # Execute adh
                print('Running AdH...')

                # Open the log file
                o_adhLog = open('adh.out', 'w')

                # Call the process
                o_adhProcess = subprocess.Popen(['mpiexec', '-n', str(self.numberOfProcessors), s_executableFullPathAdh,
                                                 filename], shell=False, stderr=subprocess.PIPE, stdout=o_adhLog)

                # Flush and close the log file
                o_adhLog.flush()
                o_adhLog.close()

                # Get the pre-AdH errors
                o_error = o_adhProcess.communicate()

                # Request the status code from the process
                i_status = o_adhProcess.returncode

                # Throw an error if the simulation did not complete successfully
                if i_status != 0:
                    # Error status received.
                    raise IOError('AdH failed to solve ' + filename + ' simulation.')

            else:
                # Solve on a single core
                # Set executable paths based on machine architecture
                if s_operatingSystem == 'Windows' and architecture == '64bit':
                    # Construct the path to the Windows 64bit executables
                    if self.version == 4.5:
                        s_executableFullPathPre = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                               'pre_adh_V4.5-WIN64.exe')
                        s_executableFullPathAdh = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                               'adh_V4.5-WIN64.exe')
                    elif self.version == 5.0:
                        s_executableFullPathPre = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                               'pre_adh_V5.0-WIN64.exe')
                        s_executableFullPathAdh = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                               'adh_V5.0-WIN64.exe')
                    else:
                        raise IOError('AdH Windows 64bit not available for this version.')

                elif s_operatingSystem == 'Windows' and architecture == '32bit':
                    # Construct the path to the Windows 32bit executables
                    if self.version == 4.5:
                        s_executableFullPathPre = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                               'pre_adh_V4.5-WIN32.exe')
                        s_executableFullPathAdh = os.path.join(s_repositoryPath, 'executables', 'windows',
                                                               'adh_V4.5-WIN32.exe')
                    else:
                        raise IOError('AdH Windows 32bit not available for this version.')

                else:
                    raise NotImplementedError('AdH solutions on non-Windows systems has yet to be implemented')

                if self.version == 4.5:
                    # Execute pre_adh
                    # Open the log file
                    o_preAdHLog = open('pre_adh.out', 'w')

                    # Call the process
                    o_preAdHProcess = subprocess.Popen(['mpiexec', '-n', '1', s_executableFullPathPre, filename],
                                                       shell=False, stderr=subprocess.PIPE, stdout=o_preAdHLog)

                    # Flush and close the log file
                    o_preAdHLog.flush()
                    o_preAdHLog.close()

                    # Get the pre-AdH errors
                    o_error = o_preAdHProcess.communicate()

                    # Request the status code from the process
                    i_status = o_preAdHProcess.returncode

                    # Throw an error if the simulation did not complete successfully
                    if i_status != 0:
                        # Error status received.
                        raise IOError('PreAdH failed to validate ' + filename + ' simulation.')

                # Execute adh
                print('Running AdH...')

                # Open the log file
                o_adhLog = open('adh.out', 'w')

                # Call the process
                o_adhProcess = subprocess.Popen(['mpiexec', '-n', '1', s_executableFullPathAdh,
                                                 filename], shell=False, stderr=subprocess.PIPE, stdout=o_adhLog)

                # Flush and close the log file
                o_adhLog.flush()
                o_adhLog.close()

                # Get the pre-AdH errors
                o_error = o_adhProcess.communicate()

                # Request the status code from the process
                i_status = o_adhProcess.returncode

                # Throw an error if the simulation did not complete successfully
                if i_status != 0:
                    # Error status received.
                    raise IOError('AdH failed to solve ' + filename + ' simulation.')

        # Otherwise solve on HPC
        else:
            raise NotImplementedError('Capability to run on HPC not yet implemented')

        # Go back to plans directory
        os.chdir(s_curdir)

    def extract_hotstart(self, result, time):
        """
        Extracts a hotstart dataset from a given result dataset and a given time
        """
        pass

    def create_hotstart(self):
        """Create hotstart from given constants or surface"""
        pass
