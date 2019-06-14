import logging
import os

import pandas as pd
import xarray as xr

from .utils import get_file_list

log = logging.getLogger('nmutils')

ROW_TYPES = ['ND', 'E2L', 'E3L', 'E3T', 'E4T', 'E6T', 'E4Q', 'E8Q', 'E9Q']


def parse_mesh_header(file_name):
    """Parse a 2dm/3dm file to extract data structure

    Args:
        file_name(str, required):
            path to 2dm/3dm file

    Returns:
        dict with card/value pairs from header
    """
    with open(file_name, 'r') as f:
        lines = f.readlines(1024)

    header = dict()
    for index, row in enumerate(lines):
        values = row.strip().split()
        if values[0] in ROW_TYPES[1:]:
            elem_size = int(values[0][1])
            extra_cols = len(values) - elem_size - 3
            break
        else:
            header[values[0]] = ' '.join(values[1:])
    # TODO have better error statements
    header['extra_cols'] = extra_cols
    header['elem_size'] = elem_size
    return header


def get_column_names(max_element_size, num_extra_columns):
    """Generate a list of column names used for Pandas to parse the data

    Args:
        max_element_size(int, required):
            the maximum element size in the mesh elements (e.g. if the mesh only contains E3T elements then
            the max_element_size is 3. If it also contains E4Q elements then it is 4).
        num_extra_columns(int, required):
            The number of extra columns after the element nodes and the material

    Returns:
        List of column names
    """
    max_cols = max_element_size + num_extra_columns + 1  # add column for element material
    data_columns = ['cmp{}'.format(i) for i in range(max_cols)]
    names = ['row_type']
    names.extend(data_columns)
    return names


def columns_from_row_type(row_type, extra_column_names=None):
    """Generate list of column names based on the row_type card at the beginning of the line

    Args:
        row_type(str, required):
            The card at the beginning of a line, should be one of `ROW_TYPES`
        extra_column_names(list<str>, optional, default=None):
            List of column names that specify what to call extra columns if they are to be preserved

    Returns:
        List of column names
    """
    extra_column_names = extra_column_names or list()

    if row_type.upper() == 'ND':
        columns = ['x', 'y', 'z']
    else:
        num_elements = int(row_type[1])
        columns = ['v{}'.format(i) for i in range(num_elements)]
        columns.append('material')

    columns.extend(extra_column_names)
    return columns


def group_to_data_array(group_name, groups, dtype=int, name=None, extra_column_names=None, attrs=None):
    """Convert a group from a pandas.Groupby object and convert it into a xarray.DataArray

    Args:
        group_name(str, required):
            name of group to extract
        groups(pandas.Groupby, required):
            Groupby object to extract group from
        dtype(type, optional, default=int):
            type to cast values to
        name(str, optional, default=None):
            name for the DataArray, if not specified then `group_name` is used
        extra_column_names(list<str>, optional, default=None):
            list of names for extra columns
        attrs(dict, optional, default=None):
            dict of attributes to attach to the DataArray

    Returns:
        xarray.DataArray of data extracted from the `group_name` group in `groups`
    """
    if group_name not in ROW_TYPES:
        log.warning(f'Skipping rows: {group_name}')
        return

    name = name or group_name

    column_names = columns_from_row_type(group_name, extra_column_names=extra_column_names)
    group = groups.get_group(group_name)
    columns_idx = group.columns[:len(column_names)]
    values = group[columns_idx].values.astype(dtype)
    dims = name + '_ids', name + '_dims'
    coords = {dims[0]: group.index, dims[1]: column_names}

    arr = xr.DataArray(data=values, dims=dims, coords=coords, name=name, attrs=attrs)
    return arr


def parse_mesh_file(file_name, max_element_size=3, extra_column_names=None, crs=None):
    """Read in a 2dm/3dm mesh file as an xarray.Dataset object

    Args:
        file_name(str, required):
            path to 2dm/3dm file to read
        max_element_size(int, optional, default=3):
            the maximum element size in the mesh elements (e.g. if the mesh only contains E3T elements then
            the max_element_size is 3. If it also contains E4Q elements then it is 4). If the first element found in
            the file is larger than `max_element_size` then that element size will be used instead. If a larger element
            is contained later in the file then an error will be raised.
        extra_column_names(list<str>, optional, default=None):
            list of names for any extra columns that should be preserved after the element and material columns are
            parsed. If there are more columns than names given extra columns will be dropped.
        crs(cartopy.CRS, optional, default=None):
            The projection of the mesh file.

    Returns:
        xarray.Dataset with a variable for the nodes and for each element type
    """
    file_name = get_file_list(file_name)[0]
    project_name = os.path.splitext(os.path.basename(file_name))[0]
    extra_column_names = extra_column_names or list()

    header = parse_mesh_header(file_name)
    skip_rows = len(header) - 2
    num_extra_cols = max(len(extra_column_names), header.pop('extra_cols', 0))
    max_element_size = max(max_element_size, header.pop('elem_size', 0))
    names = get_column_names(max_element_size, num_extra_cols)

    all_df = pd.read_csv(file_name, delim_whitespace=True, header=None, skiprows=skip_rows,
                         names=names, index_col=1)

    data_columns = names[1:]
    group_df = all_df.groupby('row_type')[data_columns]
    keys = list(group_df.groups.keys())
    keys.remove('ND')
    nodes_arr = group_to_data_array('ND', group_df, dtype=float, name='nodes',
                                    extra_column_names=extra_column_names, attrs=header)
    mesh_ds = nodes_arr.to_dataset()
    mesh_ds.attrs['project_name'] = project_name
    mesh_ds.attrs['crs'] = '' if crs is None else crs.proj4_init

    for key in keys:
        arr = group_to_data_array(key, group_df, extra_column_names=extra_column_names)
        if arr is not None:
            mesh_ds[key] = arr

    return mesh_ds
