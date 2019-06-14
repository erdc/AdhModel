from glob import glob
from collections import namedtuple
import logging

log = logging.getLogger('nmutils')


def get_file_list(glob_string=None, file_list=None):
    """Merge a glob and an explicit list of files into one list

    Args:
        glob_string(str, optional, default=None):
            string with wild card characters to expand into a list of files with glob
        file_list(list<str>, optional, default=None):
            list of file paths

    Returns:
        merged list of file paths from expanded `glob_string` and `file_list`
    """
    file_list = file_list or list()
    files = glob(glob_string)
    files.extend(file_list)
    if len(files) == 0:
        raise ValueError('No Files were found matching {} or {}.'.format(glob_string, file_list))
    return files


def get_crs(xarr):
    """Function to get a `cartopy.CRS` instance from a parsed AdH model

    Args:
        xarr(xarray.Dataset, required):
            The parsed data from an AdH model (or mesh). Should have a `crs` attribute populated with a valid proj4 string.

    Returns:
        `cartopy.CRS` instance
    """
    from cartopy import crs as ccrs

    Ccrs_Mapping = namedtuple('Ccrs_Mapping', ['ccrs_class', 'parameter_mapping'])
    PM = namedtuple('Parameter_Mapping', ['key', 'cast_func'])

    CRS_MAP = {
        'utm': Ccrs_Mapping(ccrs.UTM, {'zone': PM('zone', int), 'south': PM('southern_hemisphere', lambda x: True)}),
        'eqc': Ccrs_Mapping(ccrs.PlateCarree, {'lon_0': PM('central_longitude', float)}),
        'merc': Ccrs_Mapping(lambda: ccrs.GOOGLE_MERCATOR, {}),
    }

    proj4_init = xarr.crs
    if not proj4_init:
        log.warning('Dataset object does not have a CRS specified.')
        return
    proj4_params = _parse_proj4(proj4_init)
    mapping = CRS_MAP[proj4_params['proj']]

    kwargs = _get_mapped_parameters(mapping.parameter_mapping, proj4_params)
    crs = mapping.ccrs_class(**kwargs)

    return crs


def _parse_proj4(proj4):
    """Parses a proj4 string into a dictionary

    Args:
        proj4(str, required):
            a valid proj4 string

    Returns:
        dict of parsed proj4 parameters
    """
    components = proj4.strip('+ ').split('+')
    proj4_params = dict()
    for c in components:
        c = c.strip().split('=')
        c.append(None)
        proj4_params[c[0]] = c[1]

    proj4_params.pop('no_defs', None)
    return proj4_params


def _get_mapped_parameters(parameter_mapping, proj4_params):
    """Returns a dict of keyword arguments for a `cartopy.CRS` instance based on proj4 parameters.

    Args:
        parameter_mapping(`Parameter_Mapping`, required):
            a named tuple containing the proj4 parameters and cooresponding keyword arguments for a `cartopy.CRS` class
        proj4_params(dict, required):
            dictionary of proj4 parameters

    Returns:
        dict of keyword arguments for a `cartopy.CRS` class
    """
    kwargs = dict()
    for key, mapping in parameter_mapping.items():
        try:
            value = proj4_params[key]
        except KeyError:
            continue
        else:
            kwargs[mapping.key] = mapping.cast_func(value)

    return kwargs
