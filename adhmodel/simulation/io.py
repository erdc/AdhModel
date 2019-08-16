# -*- coding: utf-8 -*-
import os
import pandas as pd
from .material_properties import MaterialProperties
from .material_transport_properties import MaterialTransportProperties
from .time_series import TimeSeries

__all__ = ['write_hot_start_file', 'read_bc_file',
           'write_bc_file']


def write_hot_start_file(file_name, hot_start_list):
    """
    Writes the ADH hot start file from a list of hot start data sets

    Args:
        file_name: file name for *.hot file
        hot_start_list: list of HotStartDataSet classes
    """
    with open(file_name, 'w') as mesh_file:
        for ht in hot_start_list:
            mesh_file.write('DATASET\nOBJTYPE "mesh2d"\n')
            if len(ht.values.columns) > 1:
                mesh_file.write('BEGVEC\n')
            else:
                mesh_file.write('BEGSCL\n')
            mesh_file.write('ND {}\n'.format(len(ht.values)))
            mesh_file.write('NC {}\n'.format(ht.number_of_cells))
            mesh_file.write('NAME "{}"\n'.format(ht.name))
            mesh_file.write('TS 0 0\n')
            mesh_file.write(ht.values.to_csv(sep=' ', index=False, header=False).replace('\r\n', '\n'))
            mesh_file.write('ENDDS\n')


def read_bc_file(file_name, bc_class):
    """
    Reads the *.bc file and fills the AdhModel class

    Args:
        file_name: File name of the *.bc file
        bc_class: Boundary Condition class
    """
    # set all not required to deactivated
    bc_class.operation_parameters.set_not_required(False)
    bc_class.constituent_properties.set_not_required(False)
    bc_class.model_constants.set_not_required(False)

    bc_string_cards = {'NDS', 'EGS', 'MDS', 'MTS'}
    bc_cards = {'NB', 'DB', 'BR', 'OB', 'OFF', 'WER', 'WRS', 'FLP', 'FGT', 'SLUICE', 'SLS'}
    xy_series_cards = {'XY1', 'XY2', 'XYC', 'SERIES'}
    pc_cards = {'PC', 'OC', 'OS', 'FLX', 'SOUT', 'FOUT'}
    temp_data = {}
    xy_data_list = []
    with open(file_name, "r") as file:
        for line_number, line in enumerate(file):
            # remove new line character
            line = line.rstrip()
            line_split = line.split()
            # remove blank strings
            line_split[:] = (part for part in line_split if part != '')
            # skip blank line, comment line
            if len(line_split) == 0 or line_split[0] == '' or line_split[0][0] == '!':
                continue

            try:
                if line_split[0] == 'OP':
                    read_op_cards(line_split, bc_class, temp_data)
                elif line_split[0] == 'IP':
                    read_ip_cards(line_split, bc_class, temp_data)
                elif line_split[0] == 'CN':
                    read_cn_cards(line_split, bc_class, temp_data)
                elif line_split[0] == 'MP':
                    read_mp_cards(line_split, bc_class)
                elif line_split[0] in bc_string_cards:
                    read_bc_string_cards(line_split, temp_data)
                elif line_split[0] in xy_series_cards:
                    read_xy_cards(line_split, temp_data)
                elif line_split[0] == 'FR':
                    read_fr_cards(line_split, temp_data)
                elif line_split[0] in pc_cards:
                    read_pc_cards(line_split, bc_class, temp_data)
                elif line_split[0] in bc_cards:
                    read_bc_cards(line_split, bc_class, temp_data)
                elif line_split[0] == 'TC':
                    read_tc_cards(line_split, bc_class)
                elif 'xy_type' in temp_data:
                    xyt = temp_data['xy_type']
                    if xyt == 'SERIES AWRITE':
                        labels = ['START_TIME', 'END_TIME', 'TIME_STEP_SIZE', 'UNITS']
                        xy_data_list.append([float(line_split[0]), float(line_split[1]), float(line_split[2]),
                                             int(line_split[3])])
                    elif xyt == 'SERIES WIND' or xyt == 'SERIES WAVE':
                        labels = ['X', 'Y', 'Y2']
                        xy_data_list.append([float(line_split[0]), float(line_split[1]), float(line_split[2])])
                    else:
                        labels = ['X', 'Y']
                        xy_data_list.append([float(line_split[0]), float(line_split[1])])

                    # set the time step option in the output control if we read 'SERIES DT'
                    if xyt == 'SERIES DT':
                        bc_class.time_control.time_step_option = 'Time step series (SERIES DT)'
                        bc_class.time_control.max_time_step_size_time_series = temp_data['xy_id']
                    if len(xy_data_list) == temp_data['xy_number_points']:
                        ts = TimeSeries()
                        ts.series_type = xyt
                        if xyt == 'SERIES AWRITE':
                            # objs = list(bc_class.output_control.param.output_control_option.get_range())
                            bc_class.output_control.output_control_option = 'Specify autobuild (SERIES AWRITE)'
                        ts.units = temp_data['xy_units']
                        ts.output_units = temp_data['xy_output_units']
                        ts.time_series = pd.DataFrame.from_records(xy_data_list, columns=labels)
                        if 'xy_x_location' in temp_data:
                            ts.x_location = temp_data['xy_x_location']
                            ts.y_location = temp_data['xy_y_location']
                            temp_data.pop('xy_x_location')
                            temp_data.pop('xy_y_location')
                        xy_data_list = []
                        # set time series ID as both the key and in the ID column
                        ts.series_id = temp_data['xy_id']
                        bc_class.time_series[temp_data['xy_id']] = ts
                        # empty out temp_data  #todo poor practice
                        temp_data.pop('xy_number_points')
                        temp_data.pop('xy_id')
                        temp_data.pop('xy_type')
                        temp_data.pop('xy_units')
                        temp_data.pop('xy_output_units')
            except:
                msg = 'Error reading line {} of file: {}.\nLine: {}'.format(line_number+1,
                                                                              os.path.basename(file_name), line)
                raise IOError(msg)

    lists_to_data_frames(bc_class, temp_data)


def lists_to_data_frames(bc_class, temp_data):
    """
    Converts temporary lists to DataFrames in the AdhModel class
    Args:
        bc_class: The ADH boundary condition class that holds the data
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    if 'bc_string_list' in temp_data:
        labels = ['CARD', 'ID', 'ID_0', 'ID_1']
        df = pd.DataFrame.from_records(temp_data['bc_string_list'], columns=labels)
        for x in range(1, len(labels)):
            df[labels[x]] = df[labels[x]].astype(dtype='Int64')
        bc_class.boundary_strings = df
    if 'bc_list' in temp_data:
        labels = ['CARD', 'CARD_2', 'STRING_ID', 'XY_ID1', 'XY_ID2', 'XY_ID3']
        df = pd.DataFrame.from_records(temp_data['bc_list'], columns=labels)
        for x in range(2, len(labels)):
            df[labels[x]] = df[labels[x]].astype(dtype='Int64')
        bc_class.solution_controls = df
    if 'nb_sdr_list' in temp_data:
        labels = ['CARD', 'CARD_1', 'S_ID', 'COEF_A', 'COEF_B', 'COEF_C', 'COEF_D', 'COEF_E']
        df = pd.DataFrame.from_records(temp_data['nb_sdr_list'], columns=labels)
        bc_class.stage_discharge_boundary = df
    if 'fr_list' in temp_data:
        labels = ['CARD', 'CARD_2', 'STRING_ID', 'REAL_01', 'REAL_02', 'REAL_03', 'REAL_04', 'REAL_05']
        df = pd.DataFrame.from_records(temp_data['fr_list'], columns=labels)
        bc_class.friction_controls = df
    if 'br_list' in temp_data:
        labels = ['CARD', 'CARD_1', 'C_0', 'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8']
        df = pd.DataFrame.from_records(temp_data['br_list'], columns=labels)
        bc_class.breach_controls = df
    if 'wrs_list' in temp_data:
        labels = ['CARD', 'NUMBER', 'S_UPSTREAM', 'S_DOWNSTREAM', 'WS_UPSTREAM', 'WS_DOWNSTREAM', 'LENGTH',
                  'CREST_ELEV', 'HEIGHT']
        df = pd.DataFrame.from_records(temp_data['wrs_list'], columns=labels)
        bc_class.weirs = df
    if 'fgt_list' in temp_data:
        labels = ['CARD', 'NUMBER', 'USER', 'S_UPSTREAM', 'S_DOWNSTREAM', 'FS_UPSTREAM', 'FS_DOWNSTREAM', 'COEF_A',
                  'COEF_B', 'COEF_C', 'COEF_D', 'COEF_E', 'COEF_F', 'LENGTH']
        df = pd.DataFrame.from_records(temp_data['fgt_list'], columns=labels)
        bc_class.flap_gates = df
    if 'sls_list' in temp_data:
        labels = ['CARD', 'NUMBER', 'S_UPSTREAM', 'S_DOWNSTREAM', 'SS_UPSTREAM', 'SS_DOWNSTREAM', 'LENGTH',
                  'TS_OPENING']
        df = pd.DataFrame.from_records(temp_data['sls_list'], columns=labels)
        bc_class.sluice_gates = df
    if 'cn_con_list' in temp_data:
        labels = ['CARD', 'CARD_1', 'ID', 'CONC']
        df = pd.DataFrame.from_records(temp_data['cn_con_list'], columns=labels)
        bc_class.constituent_properties.general_constituents = df
    if 'cn_snd_list' in temp_data:
        labels = ['CARD', 'CARD_1', 'ID', 'C_0', 'C_1', 'C_2', 'C_3']
        df = pd.DataFrame.from_records(temp_data['cn_snd_list'], columns=labels)
        bc_class.constituent_properties.sediment = df
    if 'flx_list' in temp_data:
        labels = ['CARD', 'S_ID']
        df = pd.DataFrame.from_records(temp_data['flx_list'], columns=labels)
        bc_class.output_control.output_flow_strings = df


def read_op_cards(line_split, bc_class, temp_data):
    """
    Reads the OP cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    try:
        card = line_split[1]
        if card == 'SW2':
            bc_class.operation_parameters.physics = 'SW2'
        elif card == 'SW3':
            bc_class.operation_parameters.physics = 'SW3'
        elif card == 'INC':
            bc_class.operation_parameters.incremental_memory = int(line_split[2])
        elif card == 'TRN':
            bc_class.operation_parameters.transport = int(line_split[2])
        elif card == 'BLK':
            bc_class.operation_parameters.blocks_per_processor = int(line_split[2])
        elif card == 'PRE':
            bc_class.operation_parameters.preconditioner_type = int(line_split[2])
        elif card == 'BT':
            bc_class.operation_parameters.vessel = True
        elif card == 'BTS':
            bc_class.operation_parameters.vessel_entrainment = True
        elif card == 'TEM':
            bc_class.operation_parameters.second_order_temporal_coefficient_active = True
            bc_class.operation_parameters.second_order_temporal_coefficient = float(line_split[2])
        elif card == 'TPG':
            bc_class.operation_parameters.petrov_galerkin_coefficient_active = True
            bc_class.operation_parameters.petrov_galerkin_coefficient = float(line_split[2])
        elif card == 'NF2':
            bc_class.operation_parameters.velocity_gradient = True
        elif card == 'WND':
            bc_class.operation_parameters.wind = True
        elif card == 'WAV':
            bc_class.operation_parameters.wave = True
        elif card == 'DAM':
            bc_class.operation_parameters.dam = True
        elif card == 'DIF':
            bc_class.operation_parameters.diffusive_wave = True
    except:
        raise IOError("Error reading OP card from *.bc file.")


def read_ip_cards(line_split, bc_class, temp_data):
    """
    Reads the OP cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    try:
        card = line_split[1]
        if card == 'NIT':
            bc_class.iteration_parameters.non_linear_iterations = int(line_split[2])
        elif card == 'NTL':
            bc_class.iteration_parameters.non_linear_residual_tolerance = float(line_split[2])
            temp_data['IP NTL'] = True
            objects = list(bc_class.iteration_parameters.param.non_linear_tolerance_option.get_range())
            if 'IP ITL' in temp_data:
                bc_class.iteration_parameters.non_linear_tolerance_option = objects[0]
            else:
                bc_class.iteration_parameters.non_linear_tolerance_option = objects[1]
        elif card == 'ITL':
            bc_class.iteration_parameters.non_linear_incremental_tolerance = float(line_split[2])
            temp_data['IP ITL'] = True
            objects = list(bc_class.iteration_parameters.param.non_linear_tolerance_option.get_range())
            if 'IP NTL' in temp_data:
                bc_class.iteration_parameters.non_linear_tolerance_option = objects[0]
            else:
                bc_class.iteration_parameters.non_linear_tolerance_option = objects[2]
        elif card == 'MIT':
            bc_class.iteration_parameters.linear_iterations = int(line_split[2])
    except:
        raise IOError("Error reading IP card from *.bc file.")


def read_cn_cards(line_split, bc_class, temp_data):
    """
    Reads the CN cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    try:
        card = line_split[1]
        if card == 'CON':
            bc_class.constituent_properties.general_constituents = \
                bc_class.constituent_properties.general_constituents.append({'ID': int(line_split[2]),
                                                                              'CONC': float(line_split[3])},
                                                                             ignore_index=True)
            # bc_class.constituent_properties.general_constituents.astype({'ID': int, 'CONC': float}, copy=False)
        elif card == 'SND':
            bc_class.constituent_properties.sand = \
                bc_class.constituent_properties.sand.append({'ID': int(line_split[2]), 'C_0': float(line_split[3]),
                                                              'C_1': float(line_split[4]), 'C_2': float(line_split[5]),
                                                              'C_3': float(line_split[6])},
                                                             ignore_index=True)
            # bc_class.constituent_properties.sand.astype({'ID': int, 'C_0': float, 'C_1': float, 'C_2': float,
            #                                               'C_3': float}, copy=False)
        elif card == 'CLA':
            bc_class.constituent_properties.clay = \
                bc_class.constituent_properties.clay.append({'ID': int(line_split[2]), 'C_0': float(line_split[3]),
                                                              'C_1': float(line_split[4]), 'C_2': float(line_split[5]),
                                                              'C_3': float(line_split[6])},
                                                             ignore_index=True)
            # bc_class.constituent_properties.clay.astype({'ID': int, 'C_0': float, 'C_1': float, 'C_2': float,
            #                                               'C_3': float}, copy=False)
        elif card == 'SAL':
            bc_class.constituent_properties.salinity = True
            bc_class.constituent_properties.salinity_id = int(line_split[2])
            bc_class.constituent_properties.reference_concentration = float(line_split[3])
        elif card == 'TMP':
            bc_class.constituent_properties.temperature = True
            bc_class.constituent_properties.temperature_id = int(line_split[2])
            bc_class.constituent_properties.reference_temperature = float(line_split[3])
            val = int(line_split[4])
            if val:
                bc_class.constituent_properties.air_water_heat_transfer = True
        elif card == 'VOR':
            bc_class.constituent_properties.vorticity = True
            bc_class.constituent_properties.vorticity_id = int(line_split[2])
            bc_class.constituent_properties.vorticity_normalization = float(line_split[3])
            bc_class.constituent_properties.vorticity_as_term = float(line_split[4])
            bc_class.constituent_properties.vorticity_ds_term = float(line_split[5])

    except:
        raise IOError("Error reading CN card from *.bc file.")


def read_mp_cards(line_split, bc_class):
    """
    Reads the MP cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
    """
    try:
        card = line_split[1]
        if card == 'MU':
            bc_class.model_constants.kinematic_viscosity = float(line_split[2])
        elif card == 'G':
            bc_class.model_constants.gravity = float(line_split[2])
        elif card == 'MUC':
            bc_class.model_constants.mannings_unit_constant = float(line_split[2])
        elif card == 'RHO':
            bc_class.model_constants.density = float(line_split[2])
        elif card == 'DTL':
            bc_class.model_constants.enable_wetting_drying = True
            bc_class.model_constants.wet_dry_limit = float(line_split[2])
        else:
            # get material ID
            if card != 'WND':
                material_id = int(line_split[2])
            else:
                material_id = int(line_split[3])

            if material_id in bc_class.material_properties:
                material_property = bc_class.material_properties[material_id]
            else:
                # get material property class
                material_property = MaterialProperties()
                # set not required to deactivated
                material_property.set_not_required(False)
                bc_class.material_properties[material_id] = material_property
            if card == 'EVS':
                material_property.eddy_viscosity_method = 'Constant (EVS)'
                material_property.vxx_eddy_viscosity = float(line_split[3])
                material_property.vyy_eddy_viscosity = float(line_split[4])
                material_property.vxy_eddy_viscosity = float(line_split[5])
            elif card == 'EEV':
                material_property.eddy_viscosity_method = 'Estimated (EEV)'
                material_property.estimated_eddy_viscosity_weighting_factor = float(line_split[3])
                material_property.estimated_eddy_viscosity_method = int(line_split[4])
            elif card == 'COR':
                material_property.coriolis = True
                material_property.coriolis_latitude = float(line_split[3])
            elif card == 'ML':
                material_property.max_refinement_level = int(line_split[3])
            elif card == 'SRT':
                material_property.refinement_tolerance = float(line_split[3])
            elif card == 'DF' or card == 'TRT':
                constituent_id = int(line_split[3])
                transport_property = material_property.transport_properties.setdefault(constituent_id,
                                                                                       MaterialTransportProperties())
                if card == 'TRT':
                    transport_property.refinement_tolerance = float(line_split[4])
                else:
                    transport_property.turbulent_diffusion_rate = float(line_split[4])
            elif card == 'WND':
                if line_split[2] == 'STR':
                    material_property.wind_properties.stress_formulation = int(line_split[4])
                elif line_split[2] == 'ATT':
                    material_property.wind_properties.attenuation = float(line_split[4])
    except:
        raise IOError("Error reading MP card from *.bc file.")


def read_xy_cards(line_split, temp_data):
    """
    Reads the xy series card from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    try:
        card = line_split[0]
        # ts = TimeSeries()
        if card == 'XY1':
            temp_data['xy_type'] = 'SERIES BC'
            temp_data['xy_id'] = int(line_split[1])
            temp_data['xy_number_points'] = int(line_split[2])
            temp_data['xy_units'] = int(line_split[3])
            temp_data['xy_output_units'] = int(line_split[4])
        elif card == 'SERIES':
            temp_data['xy_type'] = '{} {}'.format(card, line_split[1])
            temp_data['xy_id'] = int(line_split[2])
            temp_data['xy_number_points'] = int(line_split[3])
            if line_split[1] == 'AWRITE':
                temp_data['xy_units'] = temp_data['xy_output_units'] = int(line_split[4])
            elif line_split[1] == 'WIND' or line_split[1] == 'WAVE':
                temp_data['xy_x_location'] = float(line_split[4])
                temp_data['xy_y_location'] = float(line_split[5])
                temp_data['xy_units'] = int(line_split[6])
                temp_data['xy_output_units'] = int(line_split[7])
            else:
                temp_data['xy_units'] = int(line_split[4])
                temp_data['xy_output_units'] = int(line_split[5])
    except:
        raise IOError("Error reading XY card from *.bc file.")


def read_bc_string_cards(line_split, temp_data):
    """
    Reads the NDS, EGS, MDS, MTS cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    bc_string_list = temp_data.setdefault('bc_string_list', [])
    record = [float('NaN') for i in range(4)]
    record[0] = line_split[0]  # card NDS, EGS, MDS, MTS
    try:
        record[1] = int(line_split[1])   # string id
        record[2] = int(line_split[2])  # node id, cell id, material id
        if line_split[0] == 'EGS' or line_split[0] == 'MDS':
            record[3] = int(line_split[3])  # node id, cell id
        bc_string_list.append(record)
    except:
        raise IOError("Error reading boundary string from *.bc file.")


def read_fr_cards(line_split, temp_data):
    """
    Reads the FR cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    fr_list = temp_data.setdefault('fr_list', [])
    try:
        num_vals_dict = {'SAV': 1, 'URV': 2, 'EDO': 4, 'ICE': 2, 'BRD': 1}
        record = [float('NaN') for i in range(8)]
        record[0] = line_split[0]
        record[1] = line_split[1]
        record[2] = int(line_split[2])
        record[3] = float(line_split[3])
        if record[1] in num_vals_dict:
            for x in range(num_vals_dict[record[1]]):
                record[4+x] = float(line_split[4+x])
    except:
        raise IOError("Error reading friction control from *.bc file.")

    fr_list.append(record)


def read_pc_cards(line_split, bc_class, temp_data):
    """
    Reads the PC cards from the *.bc file
    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
    """
    try:
        card = line_split[1]
        oc = bc_class.output_control
        if card == 'ADP':
            oc.print_adaptive_mesh = True
        elif card == 'ELM':
            oc.print_numerical_fish_surrogate = True
        elif card == 'LVL':
            option = int(line_split[2])
            if option == 0:
                oc.screen_output_residual = True
            elif option == 1:
                oc.screen_output_all = True
        elif card == 'MEO':
            oc.screen_output_mass_error = True
        elif line_split[0] == 'OC':
            oc.oc_time_series_id = int(line_split[1])
            objects = list(oc.param.output_control_option.get_range())
            oc.output_control_option = objects[0]
        elif line_split[0] == 'FLX':
            flx_list = temp_data.setdefault('flx_list', [])
            flx_list.append(['FLX', int(line_split[1])])
        elif line_split[0] == 'SOUT':
            if line_split[1] == 'RESID':
                oc.screen_output_residual = True
            elif line_split[1] == 'ALL':
                oc.screen_output_all = True
            elif line_split[1] == 'MERROR':
                oc.screen_output_mass_error = True
            elif line_split[1] == 'NLNODE':
                oc.screen_output_worst_nonlinear_node = True
            elif line_split[1] == 'LNODE':
                oc.screen_output_worst_linear_node = True
        elif line_split[0] == 'FOUT':
            if line_split[1] == 'WIND':
                oc.file_output_wind = True
            elif line_split[1] == 'WAVE':
                oc.file_output_wave = True
            if line_split[1] == "ADAPT":
                if line_split[2] == "GRID":
                    oc.file_output_adapted_grid = True
                if line_split[2] == "SW":
                    oc.file_output_adapted_solution = True
                if line_split[2] == "CON":
                    oc.file_output_adapted_transport = True
            if line_split[1] == "SED":
                oc.file_output_sediment = True


    except:
        raise IOError("Error reading CN card from *.bc file.")


def read_bc_cards(line_split, bc_class, temp_data):
    """
    Reads the solution controls from the *.bc file

    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
   """

    # these items are NOT stored in the AdhModel.solution_controls DataFrame
    try:
        cards = {'BR', 'WER', 'WRS', 'FLP', 'FGT', 'SLUICE', 'SLS'}
        if line_split[0] in cards:
            if line_split[0] == 'BR':
                br_list = temp_data.setdefault('br_list', [])
                if len(line_split) < 11:
                    line_split.extend([float('NaN') for i in range(len(line_split), 11)])
                br_list.append(line_split[:11])
            else:
                read_wrs_fgt_sls_cards(line_split, temp_data)
            return
        if line_split[0] == 'NB' and line_split[1] == 'SDR':
            nb_sdr = temp_data.setdefault('nb_sdr_list', [])
            nb_sdr.append(['NB', 'SDR', int(line_split[2]), float(line_split[3]), float(line_split[4]),
                           float(line_split[5]), float(line_split[6]), float(line_split[7])])
            return
        if line_split[0] == 'DB' and line_split[1] == 'RAD':
            bc_class.constituent_properties.short_wave_radiation_series = int(line_split[2])
            bc_class.constituent_properties.dew_point_temperature_series = int(line_split[3])
            return

        num_vals_dict = {'NB TRN': 1, 'NB OUT': 1, 'DB OVL': 1, 'DB OVH': 2, 'DB TRN': 1}
        bc_list = temp_data.setdefault('bc_list', [])
        record = [float('NaN') for i in range(6)]
        record[0] = line_split[0]
        if record[0] == 'OFF':
            record[1] = int(line_split[1])
        elif record[0] == 'OB':
            record[1] = line_split[1]
            record[2] = int(line_split[2])
        else:
            record[1] = line_split[1]
            card = '{} {}'.format(record[0], record[1])
            record[2] = int(line_split[2])
            record[3] = int(line_split[3])
            if card in num_vals_dict:
                for x in range(num_vals_dict[card]):
                    record[4+x] = int(line_split[4+x])
    except:
        raise IOError("Error reading solution control from *.bc file.")

    bc_list.append(record)


def read_wrs_fgt_sls_cards(line_split, temp_data):
    """
    Reads the WRS, FGT, and SLS cards from the *.bc file

    Args:
        line_split: list of strings from a parsed line of a *.bc file
        temp_data: Dictionary of data that is not stored in the ADH simulation but is needed while reading the file
   """
    try:
        if line_split[0] == 'WRS':
            wrs_list = temp_data.setdefault('wrs_list', [])
            wrs_list.append(['WRS', int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[4]),
                             int(line_split[5]), float(line_split[6]), float(line_split[7]), float(line_split[8])])
        elif line_split[0] == 'FGT':
            fgt_list = temp_data.setdefault('fgt_list', [])
            fgt_list.append(['FGT', int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[4]),
                             int(line_split[5]), int(line_split[6]), float(line_split[7]), float(line_split[8]),
                             float(line_split[9]), float(line_split[10]), float(line_split[11]), float(line_split[12]),
                             float(line_split[13])])
        if line_split[0] == 'SLS':
            sls_list = temp_data.setdefault('sls_list', [])
            sls_list.append(['SLS', int(line_split[1]), int(line_split[2]), int(line_split[3]), int(line_split[4]),
                             int(line_split[5]), float(line_split[6]), int(line_split[7])])
    except:
        raise IOError("Error reading {} card from *.bc file.".format(line_split[0]))


def read_tc_cards(line_split, bc_class):
    """
    Reads the TC cards from the *.bc file

    Args:
        line_split: list of strings from a parsed line of a *.bc file
        bc_class: The ADH simulation class that will hold this information
    """
    try:
        card = line_split[1]
        if card == 'T0':
            bc_class.time_control.start_time = float(line_split[2])
            if len(line_split) > 3:
                option = int(line_split[3])
                objects = list(bc_class.time_control.param.start_time_units.get_range())
                bc_class.time_control.start_time_units = objects[option]
        elif card == 'IDT':
            objects = list(bc_class.time_control.param.time_step_option.get_range())
            bc_class.time_control.time_step_option = objects[1]
            ts_id = int(line_split[2])
            bc_class.time_control.max_time_step_size_time_series = ts_id
            # if this time series exists then set the type to SERIES DT
            if ts_id in bc_class.time_series:
                bc_class.time_series[ts_id].series_type = 'SERIES DT'
        elif card == 'TF':
            bc_class.time_control.end_time = float(line_split[2])
            if len(line_split) > 3:
                option = int(line_split[3])
                objects = list(bc_class.time_control.param.end_time_units.get_range())
                bc_class.time_control.end_time_units = objects[option]
        elif card == 'ATF':
            objects = list(bc_class.time_control.param.time_step_option.get_range())
            bc_class.time_control.time_step_option = objects[2]
            bc_class.time_control.auto_time_step_find_min_time_step_size = float(line_split[2])
            bc_class.time_control.auto_time_step_find_max_time_step_size_series = int(line_split[3])
        elif card == 'STD':
            objects = list(bc_class.time_control.param.time_step_option.get_range())
            bc_class.time_control.time_step_option = objects[0]
            bc_class.time_control.steady_state_min_time_step_size = float(line_split[2])
            bc_class.time_control.steady_state_max_time_step_size = float(line_split[3])
    except Exception as ex:
        raise IOError("Error reading TC card from *.bc file.")


def write_bc_file(file_name, bc_class, validate=False):
    """
    Writes a *.bc simulation give the information in the AdhModel class

    Args:
        file_name: File name
        bc_class: The ADH simulation class that holds all simulation information
        validate: Run validation on parameters
    """
    if validate:
        if bc_class.operation_parameters.physics == 'SW2':
            bc_class.validate_sw2()
        else:
            raise RuntimeWarning('Validation only available for SW2')

    with open(file_name, "w") as bc_file:
        bc_file.write('! ADH BC File - written by adhparam\n\n')
        write_op_cards(bc_file, bc_class)
        write_ip_cards(bc_file, bc_class)
        write_cn_cards(bc_file, bc_class)
        write_mp_cards(bc_file, bc_class)
        write_boundary_string_cards(bc_file, bc_class)
        write_time_series_cards(bc_file, bc_class)
        write_fr_cards(bc_file, bc_class)
        write_pc_cards(bc_file, bc_class)
        write_solution_control_cards(bc_file, bc_class)
        write_tc_cards(bc_file, bc_class)
        bc_file.write('END\n')


def number_of_constituents(bc_class):
    """
    Calculates the number of constituents

    Args:
        bc_class: The ADH simulation class that holds all simulation information

    Returns:
        The number of transport constituents

    """
    num_trn = 0
    cn = bc_class.constituent_properties
    if cn.salinity:
        num_trn += 1
    if cn.temperature:
        num_trn += 1
    if cn.vorticity:
        num_trn += 1
    if not cn.general_constituents.empty:
        num_trn += len(cn.general_constituents.index)
    if not cn.sand.empty:
        num_trn += len(cn.sand.index)
    if not cn.clay.empty:
        num_trn += len(cn.clay.index)
    return num_trn


def write_op_cards(bc_file, bc_class):
    """
    Writes the OP cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    op = bc_class.operation_parameters
    bc_file.write('! Operation Parameters\n')
    bc_file.write('OP {}\n'.format(op.physics))
    bc_file.write('OP INC {}\n'.format(op.incremental_memory))
    bc_file.write('OP TRN {}\n'.format(op.transport))
    bc_file.write('OP BLK {}\n'.format(op.blocks_per_processor))
    bc_file.write('OP PRE {}\n'.format(op.preconditioner_type))
    if op.vessel:
        bc_file.write('OP BT\n')
        if op.vessel_entrainment:
            bc_file.write('OP BTS\n')
    if op.second_order_temporal_coefficient_active:
        bc_file.write('OP TEM {}\n'.format(op.second_order_temporal_coefficient))
    if op.petrov_galerkin_coefficient_active:
        bc_file.write('OP TPG {}\n'.format(op.petrov_galerkin_coefficient))
    if op.velocity_gradient:
        bc_file.write('OP NF2\n')
    if op.wind:
        bc_file.write('OP WND\n')
    if op.wave:
        bc_file.write('OP WAV\n')
    if op.dam:
        bc_file.write('OP DAM\n')
    if op.diffusive_wave:
        bc_file.write('OP DIF\n')

    bc_file.write('\n') # blank line at the end of the Operation Parameters


def write_ip_cards(bc_file, bc_class):
    """
    Writes the IP cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    ip = bc_class.iteration_parameters
    bc_file.write('! Iteration Parameters\n')
    bc_file.write('IP NIT {}\n'.format(ip.non_linear_iterations))
    index = list(ip.param.non_linear_tolerance_option.get_range()).index(ip.non_linear_tolerance_option)
    if index == 0 or index == 1:
        bc_file.write('IP NTL {}\n'.format(ip.non_linear_residual_tolerance))
    if index == 0 or index == 2:
        bc_file.write('IP ITL {}\n'.format(ip.non_linear_incremental_tolerance))
    bc_file.write('IP MIT {}\n'.format(ip.linear_iterations))

    bc_file.write('\n') # blank line at the end of the Iteration Parameters


def write_cn_cards(bc_file, bc_class):
    """
    Writes the CN cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    cn = bc_class.constituent_properties
    bc_file.write('! Constituent Properties\n')
    if not cn.general_constituents.empty:
        # bc_file.write(cn.general_constituents.to_csv(sep=' ', index=False, header=False).replace('\r\n', '\n'))
        for index, row in bc_class.constituent_properties.general_constituents.iterrows():
            bc_file.write(
                'CN CON {} {}\n'.format(row['ID'].astype('int'), row['CONC']))
    if not cn.sand.empty:
        # bc_file.write(cn.sand.to_csv(sep=' ', index=False, header=False).replace('\r\n', '\n'))
        for index, row in bc_class.constituent_properties.sand.iterrows():
            bc_file.write(
                'CN SND {} {} {} {} {}\n'.format(row['ID'].astype('int'), *row[['C_0', 'C_1', 'C_2', 'C_3']].values))
    if not cn.clay.empty:
        # bc_file.write(cn.clay.to_csv(sep=' ', index=False, header=False).replace('\r\n', '\n'))
        for index, row in bc_class.constituent_properties.clay.iterrows():
            bc_file.write(
                'CN CLA {} {} {} {} {}\n'.format(row['ID'].astype('int'), *row[['C_0', 'C_1', 'C_2', 'C_3']].values))
    if cn.salinity:
        bc_file.write('CN SAL {} {}\n'.format(cn.salinity_id, cn.reference_concentration))
    if cn.temperature:
        bc_file.write('CN TMP {} {}\n'.format(cn.temperature_id, cn.reference_temperature))
    if cn.vorticity:
        bc_file.write('CN VOR {} {} {} {}\n'.format(cn.vorticity_id, cn.vorticity_normalization,
                                                    cn.vorticity_as_term, cn.vorticity_ds_term))

    bc_file.write('\n')  # blank line at the end of the Constituent Properties


def write_mp_cards(bc_file, bc_class):
    """
    Writes the MP cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    mc = bc_class.model_constants
    bc_file.write('! Global Material Properties\n')
    bc_file.write('MP MU {}\n'.format(mc.kinematic_viscosity))
    bc_file.write('MP G {}\n'.format(mc.gravity))
    bc_file.write('MP RHO {}\n'.format(mc.density))
    if mc.enable_wetting_drying:
        bc_file.write('MP DTL {}\n'.format(mc.wet_dry_limit))
    bc_file.write('MP MUC {}\n'.format(mc.mannings_unit_constant))
    bc_file.write('\n')

    bc_file.write('! Material Properties\n')
    for id, mat_prop in bc_class.material_properties.items():
        if mat_prop.eddy_viscosity_method == 'Constant (EVS)':
            bc_file.write('MP EVS {} {} {} {}\n'.format(id, mat_prop.vxx_eddy_viscosity, mat_prop.vyy_eddy_viscosity,
                                                        mat_prop.vxy_eddy_viscosity))

        elif mat_prop.eddy_viscosity_method == 'Estimated (EEV)':
            objects = list(mat_prop.param.estimated_eddy_viscosity_method.get_range())
            bc_file.write('MP EEV {} {} {}\n'.format(id, mat_prop.estimated_eddy_viscosity_weighting_factor,
                                                     mat_prop.estimated_eddy_viscosity_method))
        if mat_prop.coriolis:
            bc_file.write('MP COR {} {}\n'.format(id, mat_prop.coriolis_latitude))
        bc_file.write('MP SRT {} {}\n'.format(id, mat_prop.refinement_tolerance))
        bc_file.write('MP ML {} {}\n'.format(id, mat_prop.max_refinement_level))
        if bc_class.operation_parameters.transport != 0:
            for id1, tran_prop in mat_prop.transport_properties.items():
                bc_file.write('MP TRT {} {} {}\n'.format(id, id1, tran_prop.refinement_tolerance))
                if tran_prop.refinement_tolerance > 0 and mat_prop.eddy_viscosity_method == 'Constant (EVS)':
                    bc_file.write('MP DF {} {} {}\n'.format(id, id1, tran_prop.turbulent_diffusion_rate))
        if bc_class.operation_parameters.wind:
            wnd = mat_prop.wind_properties
            bc_file.write('MP WND STR {} {}\n'.format(id, wnd.stress_formulation))
            bc_file.write('MP WND ATT {} {}\n'.format(id, wnd.attenuation))

    bc_file.write('\n')  # blank line at the end of the Material Properties


def write_boundary_string_cards(bc_file, bc_class):
    """
    Writes the NDS, EGS, MDS, MTS cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    bs = bc_class.boundary_strings
    if not bs.empty:
        bc_file.write('! Boundary Strings\n')
        bc_file.write(bs.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))
        bc_file.write('\n')  # blank line after Boundary Strings


def write_time_series(bc_file, ts, series_id):
    """
    Writes a time series to the *.bc file

    Args:
        bc_file: the *.bc file
        ts: TimeSeries class that is written to the file
        series_id: ID of the time series
    """
    # construct the card line with series type, series number, and number of entries
    line = '{} {} {}'.format(ts.series_type, series_id, len(ts.time_series.index))
    # for output series, add the output units
    if ts.series_type == 'SERIES AWRITE':
        line = '{} {}'.format(line, ts.output_units)  # todo this looks wrong
    # for wind and wave series, add the location and input/output units
    elif ts.series_type == 'SERIES WIND' or ts.series_type == 'SERIES WAVE':
        line = '{} {} {} {} {}'.format(line, ts.x_location, ts.y_location, ts.units,
                                       ts.output_units)
    else:
        # for all other series types, add input/output units
        line = '{} {} {}'.format(line, ts.units, ts.output_units)
    # write the constructed line
    bc_file.write('{}\n'.format(line))
    bc_file.write(ts.time_series.to_csv(sep=' ', index=False, header=False, ).replace('\r\n', '\n'))
    bc_file.write('\n')


def write_time_series_cards(bc_file, bc_class):
    """
    Writes the time series cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    # add header for the time series section
    bc_file.write('! Time Series\n')
    awrite_key = -1
    dt_key = -1

    for key, ts in bc_class.time_series.items():
        if ts.series_type == 'SERIES AWRITE':
            # store the output series number
            awrite_key = key
        elif ts.series_type == 'SERIES DT':
            # store the timestep series number
            dt_key = key
        else:
            # write all other series
            write_time_series(bc_file, ts, key)
    bc_file.write('\n')  # blank line after Time Series

    # write the time step series
    if dt_key != -1:
        # write header for time step series section
        bc_file.write('! Time step time series\n')
        write_time_series(bc_file, bc_class.time_series[dt_key], dt_key)
        bc_file.write('\n')  # blank line after Time step time series

    # write the output series
    if awrite_key != -1:
        # write header for time step series
        bc_file.write('! Output series\n')
        write_time_series(bc_file, bc_class.time_series[awrite_key], awrite_key)
        bc_file.write('\n')  # blank line after Output series


def write_fr_cards(bc_file, bc_class):
    """
    Writes the FR cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    fr = bc_class.friction_controls
    if not fr.empty:
        bc_file.write('! Friction Controls\n')
        bc_file.write(fr.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))
        bc_file.write('\n')  # blank line after Friction Controls


def write_pc_cards(bc_file, bc_class):
    """
    Writes the PC cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    bc_file.write('! Output Control\n')
    oc = bc_class.output_control
    objects = list(oc.param.output_control_option.get_range())
    if oc.output_control_option == objects[0]:
        bc_file.write('OC {}\n'.format(oc.oc_time_series_id))
    ofs = oc.output_flow_strings
    if not ofs.empty:
        bc_file.write(ofs.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))

    if oc.print_adaptive_mesh:
        bc_file.write('PC ADP\n')
    if oc.print_numerical_fish_surrogate:
        bc_file.write('PC ELM\n')
    if oc.screen_output_residual:
        bc_file.write('SOUT RESID\n')
    if oc.screen_output_all:
        bc_file.write('SOUT ALL\n')
    if oc.screen_output_mass_error:
        bc_file.write('SOUT MERROR\n')
    if oc.screen_output_worst_nonlinear_node:
        bc_file.write('SOUT NLNODE\n')
    if oc.screen_output_worst_linear_node:
        bc_file.write('SOUT LNODE\n')
    if oc.file_output_wind:
        bc_file.write('FOUT WIND\n')
    if oc.file_output_wave:
        bc_file.write('FOUT WAVE\n')
    if oc.file_output_adapted_grid:
        bc_file.write('FOUT ADAPT GRID\n')
    if oc.file_output_adapted_solution:
        bc_file.write('FOUT ADAPT SW\n')
    if oc.file_output_adapted_transport:
        bc_file.write('FOUT ADAPT CON\n')
    if oc.file_output_sediment:
        bc_file.write('FOUT SED\n')

    bc_file.write('\n')  # blank line after Output Control


def write_solution_control_cards(bc_file, bc_class):
    """
    Writes the NB, DB, BR, WER, FLP, SLUICE cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    bc_file.write('! Solution Controls\n')
    sc = bc_class.solution_controls
    if not sc.empty:
        bc_file.write(sc.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))

    nb_sdr = bc_class.stage_discharge_boundary
    if not nb_sdr.empty:
        bc_file.write(nb_sdr.to_csv(sep=' ', na_rep='', index=False, header=False, ).replace('\r\n', '\n'))

    if bc_class.constituent_properties.temperature:
        cp = bc_class.constituent_properties
        bc_file.write('DB RAD {} {}\n'.format(cp.short_wave_radiation_series, cp.dew_point_temperature_series))

    write_breach_control_cards(bc_file, bc_class)

    wer = bc_class.weirs
    if not wer.empty:
        bc_file.write('WER {}\n'.format(len(wer.index)))
        bc_file.write(wer.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))

    flp = bc_class.flap_gates
    if not flp.empty:
        bc_file.write('FLP {}\n'.format(len(flp.index)))
        bc_file.write(flp.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))

    sls = bc_class.sluice_gates
    if not sls.empty:
        bc_file.write('SLUICE {}\n'.format(len(sls.index)))
        bc_file.write(sls.to_csv(sep=' ', na_rep='', index=False, header=False,).replace('\r\n', '\n'))

    bc_file.write('\n')  # blank line after Solution Controls


def write_breach_control_cards(bc_file, bc_class):
    """
    Writes the BR cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    br = bc_class.breach_controls
    if br.empty:
        return
    format_dict = {'JAI': [str, str, int,     int, float, float, float, float,   int,   int],
                   'SAS': [str, str, int,     int, float, float, float, float,   int,   int],
                   'MLM': [str, str, float, float, float, float],
                   'FRO': [str, str, int,     int, float, float, float, float, float, int, int],
                   'BRC': [str, str, int,     int, float, float, float, float,   int, int],
                   'VTG': [str, str, int,     int, float, float, float,   int, float, int, int],
                   'FER': [str, str, int,     int, float, float, float,   int, float, int, int],
                   'USR': [str, str, int,     int]
                   }
    try:
        col_list = ['CARD', 'CARD_1', 'C_0', 'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6', 'C_7', 'C_8']
        for idx, row in br.iterrows():
            card = row['CARD_1']
            format_list = format_dict[card]
            for index, item in enumerate(format_list):
                record = row[col_list[index]]
                mystr = '{} '.format(item(record))
                bc_file.write(mystr)
            bc_file.write('\n')

    except:
        raise IOError('Error writing BR card.')


def write_tc_cards(bc_file, bc_class):
    """
    Writes the TC cards to the *.bc file

    Args:
        bc_file: the *.bc file
        bc_class: The ADH simulation class that holds all simulation information
    """
    bc_file.write('! Time Controls\n')
    tc = bc_class.time_control
    objects = list(tc.param.start_time_units.get_range())
    bc_file.write('TC T0 {} {}\n'.format(tc.start_time, objects.index(tc.start_time_units)))
    objects = list(tc.param.end_time_units.get_range())
    bc_file.write('TC TF {} {}\n'.format(tc.end_time, objects.index(tc.end_time_units)))
    objects = list(tc.param.time_step_option.get_range())
    if tc.time_step_option == objects[0]:
        bc_file.write('TC STD {} {}\n'.format(tc.steady_state_min_time_step_size, tc.steady_state_max_time_step_size))
    elif tc.time_step_option == objects[2]:
        bc_file.write('TC ATF {} {}\n'.format(tc.auto_time_step_find_min_time_step_size,
                                              tc.auto_time_step_find_max_time_step_size_series))
    # TC IDT replaced by SERIES DT
    # else:
    #    bc_file.write('TC IDT {}\n'.format(tc.max_time_step_size_time_series))

    bc_file.write('\n')  # blank line after Time Controls





