import logging
import os
import re
from io import StringIO
from functools import partial

import pandas as pd
import xarray as xr

from .utils import get_file_list

log = logging.getLogger('nmutils')

COLUMNS_IN_TIMESTEP_ROW = 3
TIMESTEP_ROW_REGEX = 'TS\s+0\s+([\d\.+e]*)\n'


def get_variable_from_file_name(file_name, project_name=None):
    """Parse the name of a .dat as a substitute for the variable name

    Args:
        file_name(str, required):
            path to .dat file
        project_name(str, optional, default=None):
            The root name of a project

    Returns:
          String, portion of the file name that represents the variable
    Note:
        The assumption is that the `file_name` is in the form `./path/to/dat/file/{project_name}_{variable}.dat`
    """
    project_name = project_name or ''
    name, _ = os.path.splitext(os.path.basename(file_name))
    variable = name.split(project_name + '_')[-1]
    return variable


def parse_dat_header(file_name):
    """Parse the header of a .dat file and return as a dict

    Args:
        file_name(str, required):
            path to .dat file

    Returns:
        dict of card/value pairs from the header

    Note:
        The first line of data is also parsed to determine the dimensionality of the data
    """
    with open(file_name, 'r') as f:
        text = f.read(1024)
    return _parse_header(text)


def _parse_header(text, file_name=None):
    """Helper function for parsing .dat and .hot file headers

    Args:
        text(str, required):
            the content of the first chuck of the file
        file_name(str, optional, default=None):
            The name of the file being parsed (used to display information in error messages)

    Returns:
        dict of card/value pairs from the header

    Note:
        This function was split from the `parse_dat_header` function in order to parse multiple headers from a .hot file
    """
    file_name = file_name or text
    matches = re.split(TIMESTEP_ROW_REGEX, text)
    if len(matches) < 3:
        raise ValueError('The file {} is not formatted correctly. Skipping...'.format(file_name))
    header = [row.split() for row in matches[0].splitlines()]
    header = {row[0].upper(): (' '.join(row[1:]).strip('" ') if len(row) > 1 else '') for row in header}
    dimensionality = len(matches[2].splitlines()[0].strip().split())
    header['DIM'] = dimensionality
    return header


def parse_dat_file(file_name, project_name=None):
    """Parse a .dat file into a xarray DataArray

    Args:
        file_name(str, required):
            path to .dat file
        project_name(str, optional, default=None):
            The root name of a project (used to extract the variable from the `file_name`
            but only if the file header doesn't contain a value for `NAME`)

    Returns:
        xarray.DataArray object
    """
    file_name = get_file_list(file_name)[0]
    header = parse_dat_header(file_name)
    header.setdefault('NAME', get_variable_from_file_name(file_name, project_name=project_name))
    return _parse_data(file_name, header)


def _parse_data(file_name, header):
    """Helper function for parsing data from .dat and .hot files

    Args:
        file_name(str or file-like, required:
            path to .dat file or a StringIO with data from .hot file
        header(dict, required):
            dict containing information parsed from the header fo the file
        
    Returns:
        xarray.DataArray object
        
    Note:
        This function was split from the `parse_dat_file` function in order to parse multiple datasets from a .hot file
    """
    skip_rows = len(header) - 1
    variable = header.pop('NAME')
    col_names = list(range(max(COLUMNS_IN_TIMESTEP_ROW, header['DIM'])))
    df = pd.read_table(file_name, delim_whitespace=True, skiprows=skip_rows,
                       header=None, names=col_names, engine='c', dtype=object)

    # extract the timesteps from the dataframe
    times = df[df[0] == 'TS'][2].astype(float)

    # determine if there are one or more timesteps and calculate the number of nodes in each timestep
    if times.shape[0] > 1:
        values_shape = times.index[1] - times.index[0] - 1
        dims = ['times', 'nodes_ids']
        coords = {'times': times}
    else:
        values_shape = df.shape[0] - 1
        dims = ['init_time', 'nodes_ids']
        coords = None

    # strip out time step lines, drop nan columns, and convert values to floats
    vals_df = df[df.index % (values_shape + 1) > 0].dropna(axis=1).astype(float)

    # check that the last timestep is complete
    if vals_df.index[-1] - times.index[-1] < values_shape:
        log.warning('Timestep {} is incomplete; Skipping...'.format(times.iloc[-1]))
        vals_df = vals_df[vals_df.index < times.index[-1]]
        times = times[times.index < times.index[-1]]

    # adjust shape/dims for multi-dimentional datasets (e.g. velocity)
    new_shape = [times.shape[0], values_shape]
    if header['DIM'] > 1:
        new_shape.append(header['DIM'])
        dims.append('{}_dims'.format(variable))

    values = vals_df.values.reshape(new_shape)

    vals_arr = xr.DataArray(
        name=variable,
        coords=coords,
        dims=dims,
        data=values,
        attrs=header,
    )
    return vals_arr


def parse_dat_files(glob_string='', file_list=None, project_name=None):
    """Parse a list of .dat files into a xarray Dataset

    Args:
        glob_string(str, optional, default=''):
            A pathname with wildcard characters to read a list of files with the glob module
        file_list(list<str>, optional, default=None):
            A list of file paths to .dat files to read in
        project_name(str, optional, default=None):
            the root name of the project

    Returns:
        xarray.DataArray object with a variable for reach .dat file parsed

    Note:
        Either `glob_string` or `file_list` must be specified.
        If both are specified then the results of expanding `glob_string` are extended by the `file_list`
    """
    files = get_file_list(glob_string, file_list)
    arr_func = partial(parse_dat_file, project_name=project_name)

    return _aggregate_dataset(files, arr_func)


def parse_hot_file(file_name):
    """Parse a .hot file into an xarray Dataset

    Args:
        file_name(str, required):
            path to .hot file

    Returns:
        xarray.Dataset object with a variable for each dataset listed in the .hot file
    """
    file_name = get_file_list(file_name)[0]
    with open(file_name, 'r') as f:
        text = f.read()
    data = text.split('ENDDS\n')
    assert data[-1] == ''
    data = data[:-1]

    arr_func = partial(_parse_hot_dataset, file_name=file_name)
    return _aggregate_dataset(data, arr_func)


def _parse_hot_dataset(text, file_name):
    """Helper function for `parse_hot_file` to parse one dataset in the file

    Args:
        text(str, required):
            the contents of one dataset (including the header) from a .hot file
        file_name(str, required):
            path to the .hot file (used for error messages)

    Returns:
        xarray.DataArray object
    """

    header = _parse_header(text, file_name=file_name)
    text_io = StringIO(text)
    arr = _parse_data(text_io, header)
    return arr


def _aggregate_dataset(iterable, data_array_func):
    """Helper function for `parse_dat_files` and `parse_hot_file`. Used to aggregate xarray DataArrays into a Dataset

    Args:
        iterable(iterable, required):
            An iterable of file-like objects to parse into a DataArray
        data_array_func(func, required):
            A function that parses a file into a DataArray

    Returns:
        xarray.Dataset composed of variables from each of the files in `iterable`
    """
    d = xr.Dataset()
    for arg in iterable:
        try:
            arr = data_array_func(arg)
            append_array_to_dataset(d, arr)
        except Exception as e:
            log.warning('Cannot append array because of the following error\n{}'.format(e))
            continue

    return d


def append_array_to_dataset(dataset, array):
    if array.name in dataset:
        i = 1
        while f'{array.name}({i})' in d:
            i += 1
        array.name = f'{array.name}({i})'
    dataset[array.name] = array
