import os
import networkx as nx
from scipy.spatial import distance
import geoviews as gv
import cartopy.crs as ccrs
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path

import param
import logging

from genesis.mesh import Unstructured2D
from .simulation.simulation import AdhSimulation
from .simulation import parse_mesh_file
from .simulation.utils import get_crs

log = logging.getLogger('adh')

ROOTDIR = os.path.dirname(os.path.dirname(__file__))


class AdhMesh(Unstructured2D):
    """"""
    current_sim = param.ClassSelector(default=AdhSimulation(), class_=AdhSimulation, precedence=-1)

    result_time = param.ObjectSelector()

    def __init__(self, **params):
        super(AdhMesh, self).__init__(**params)

    def read(self, path, project_name='*', crs=ccrs.UTM(1), fmt='nc'):
        file_path = Path(path)
        if not file_path.is_file():
            file_path = file_path / f'{project_name}.{fmt}'
        if fmt == 'nc':
            xarr = xr.open_dataset(file_path)
        elif fmt == '2dm' or fmt == '3dm':
            mesh_path = os.path.join(file_path)
            xarr = parse_mesh_file(mesh_path, crs=crs)

        return self.from_xarray(xarr, crs=crs)

    def write(self, file_name, fmt='nc'):
        fmts = {
            'nc': self.to_netcdf,
            '2dm': self.to_3dm,
            '3dm': self.to_3dm,
        }
        if fmt not in fmts:
            raise ValueError(f'Unsupported format. Must be one of {fmts.keys()}')
        return fmts[fmt](file_name)

    def from_xarray(self, xarr, crs=ccrs.UTM(1)):
        self.name = xarr.nodes.attrs.get('MESHNAME').strip('\'').strip('"')
        self.verts = xarr.nodes.to_pandas()  # store as one-based
        self.tris = xarr.E3T.to_pandas() - 1  # store as zero-based (will be paired with mesh_points which uses zero-based indices)
        # TODO: why don't we just store the mesh as a xarray instead of pandas?
        file_crs = get_crs(xarr)

        if crs:
            if crs.proj4_params != file_crs.proj4_params:
                log.warning('Specified crs ({}) does not match file ({}). Defaulting to file crs'.format(
                    crs.proj4_params, file_crs.proj4_params))
                self.projection.set_crs(file_crs)

            else:
                self.projection.set_crs(crs)
        else:
            self.projection.set_crs(file_crs)

        self.reproject_points()
        self.tri_mesh = gv.TriMesh((self.tris[['v0', 'v1', 'v2']], self.mesh_points))

    def to_xarray(self):
        mesh_ds = xr.DataSet()
        mesh_ds['nodes'] = self.nodes.to_xarray()
        mesh_ds['E3T'] = self.elements.to_xarray()
        return mesh_ds

    def to_netcdf(self, file_name):
        mesh_ds = self.to_xarray()
        mesh_ds.to_netcdf(file_name)

    def to_3dm(self, file_name):
        """
        Writes a *.3dm file that defines a mesh.

        Args:
            file_name: file name of *.3dm file
            mesh: Mesh class with mesh information
        """
        import csv
        with open(file_name, 'w') as mesh_file:
            mesh_file.write('MESH2D\n')
            if len(self.name) > 0:
                mesh_file.write('MESHNAME "{0}"\n'.format(self.name))
        nodes = self.tris.copy() + 1
        nodes.insert(0, 'id', nodes.index)
        nodes.insert(0, 'card', 'E3T', allow_duplicates=True)
        nodes.to_csv(file_name, mode='a', sep=' ', index=False, header=False, line_terminator='',
                     float_format='%10.6f', quoting=csv.QUOTE_NONE, quotechar='',  escapechar=' ')
        cells = self.verts.copy()
        cells.insert(0, 'id', cells.index)
        cells.insert(0, 'card', 'ND', allow_duplicates=True)
        cells.to_csv(file_name, mode='a', sep=' ', index=False, header=False, line_terminator='',
                     float_format='%10.6f', quoting=csv.QUOTE_NONE, quotechar='',  escapechar=' ')

    def reproject_points(self):
        """Method to set the mesh points as gv points (projected if necessary)"""
        # reproject the points
        self.mesh_points = gv.operation.project_points(
             gv.Points(self.verts[['x', 'y', 'z']], vdims=['z'], crs=self.projection.get_crs()),
             projection=ccrs.GOOGLE_MERCATOR)

    def validate(self):
        """
        Checks for errors in the AdH 2D mesh
        """
        if self.verts.empty():
            log.error('Mesh nodes not set.')

        if self.tris.empty():
            log.error('Mesh elements not set.')


def check_polygon(poly, rtol=1e-05, atol=1e-08):
    """method to fix duplicate start and endpoints"""
    # if the first and last nodes in an array are equal
    if np.allclose(poly[0], poly[-1], rtol, atol, equal_nan=True):
        # strip out the last node and return
        return poly[:-1]
    else:
        # otherwise, return the original poly
        return poly


def xmsmesh_to_dataframe(pts, cells):
    """
    Convert mesh pts and cells to dataframe

    Args:
      pts (MultiPolyMesherIo.points): Points from a MultiPolyMesherIo
      cells (MultiPolyMesherIo.cells: Cells from a MultiPolyMesherIo

    Returns:
      pd.DataFrame: MultiPolyMesherIo points in a dataframe
      pd.DataFrame: MultiPolyMesherIo cells in a dataframe
    """
    r_pts = pd.DataFrame(pts, columns=['x', 'y', 'z'])
    r_cells = pd.DataFrame([(cells[x + 2], cells[x + 3], cells[x + 4]) for x in range(0, len(cells), 5)],
                           columns=['v0', 'v1', 'v2'])
    return r_pts, r_cells


def elements_by_materials(tris, mats):
    """
    Get list of elements given a list of material numbers
    Args:
        tris: zero indexed df of elements
        mats: list of material ints, zero indexed

    Returns:
        elements: list of elements, zero indexed

    """

    # create an empty array for the elements
    elements = np.array([])
    # loop over all the materials given
    for mat in mats:
        # get the subset of this material
        df = tris.loc[tris['material'] == mat]
        # get the indices of the returned df
        arr = df.index.values
        # add the indices to the list of elements
        elements = np.insert(elements, 0, arr)

    return elements


def elements_by_node(tris, nodes):
    """
    Get list of elements connected to a given node
    Args:
        tris: zero indexed df of elements
        nodes: list of zero indexed node number, int

    Returns:
        elements: set of elements zero indexed
    """
    # create an empty list for the elements
    elements = set()
    for node in nodes:
        # loop over all the materials given
        for vert in ['v0', 'v1', 'v2']:
            # get the subset of this material
            # df = tris.loc[tris['mat'] == mat]
            df = tris.loc[tris[vert] == node]
            # add the indices to the df to elements
            elements.update(df.index.values.astype(int))

    return elements


def nodes_by_elements(tris, elements):
    """
    Get list of nodes included in a list of elements
    Args:
        tris: zero indexed df of elements
        elements: zero indexed list of elements

    Returns:
        nodes zero indexed set of nodes
    """
    try:
        # get a subset with the requested elements
        df = tris.loc[elements]
    except:
        log.error('Element index not found in tris dataframe')
    # get the values as array
    nodes = set(df[['v0', 'v1', 'v2']].values.flatten())
    # return the zero-based node array
    return nodes


def reduce_by_box(verts, xmin, xmax, ymin, ymax):
    """Reduce the vertex dataframe based on a bounding box"""
    verts = verts.loc[verts['x'] > xmin]
    verts = verts.loc[verts['x'] < xmax]
    verts = verts.loc[verts['y'] > ymin]
    verts = verts.loc[verts['y'] < ymax]

    # check to ensure something was found
    if len(verts) == 0:
        log.error('No nodes were found in the bounding box: '
                  'xmin={}, xmax={}, ymin={}, ymax={}'.format(xmin, xmax, ymin, ymax))

    return verts


def nearest_node(verts, point):
    """get the nearest node to a given location
    NOTE: this is computationally expensive and
    it is recommended that the verts df be reduced
    as much as possible before calling this function."""
    # convert to array
    arr = np.array(list(zip(verts['x'], verts['y'])))
    # calculate the distance between the vertices in the bounding box
    # to the reservoir point
    d = distance.cdist(arr, point)
    # find the node at the minimium distance
    idx = np.argmin(d)
    # get the index of the nearest node
    node = verts.iloc[idx].name
    # get the distance
    dist = d[idx]

    return node, dist


def create_graph(verts, tris, dist=False):
    """
    Create networkx graph from verts and tris dfs
    Args:
        verts: Dataframe of vertices with columns for 'x' and 'y'
        tris: Dataframe of elements with columns for 'v0', 'v1', and 'v2'
        dist: Boolean toggle for calculating the distances of each edge
                (NOTE: this is computationally expensive for large meshes)

    Returns:
        g: networkx graph. Optional: Distances stored as 'dist' weight
    """
    # create empty graph
    g = nx.Graph()
    # add the nodes to the graph
    g.add_nodes_from(verts)

    # for each edge of element
    edges = [['v0', 'v1'], ['v1', 'v2'], ['v0', 'v2']]
    for edg in edges:
        # reorganize the edge vertices
        elist = list(zip(tris[edg[0]], tris[edg[1]]))

        # if the distance calculation is required
        if dist:
            try:
                # reorganize the vertices
                n1 = np.array((verts['x'][elist[0]], verts['y'][elist[0]]))
                n2 = np.array((verts['x'][elist[1]], verts['y'][elist[1]]))
                # calculate the distance between the vertices
                d = distance.cdist(n1, n2)

            except:
                d = np.ones(len(elist)) * 9999.9999

            # reorganize the list to include the distances
            newlist = list(zip(tris[edg[0]], tris[edg[1]], d))

            # add the edges from the list
            g.add_weighted_edges_from(newlist, weight='dist')
        else:
            # at the edge from the list
            g.add_edges_from(elist)

    # return the networkx graph
    return g


def get_path_by_loc(g, locs, verts):
    """
    Method to find the shortest path of mesh nodes between a series of user-specified locations.
    Units of the location points must be the same coordinate system as verts.
    Args:
        g: networkx graph - MUST include distances as 'dist' weights
        locs: list of x,y pairs of the locations
        verts: dataframe of vertices with columns 'x' and 'y'

    Returns:
        path: list of node numbers along the path
    """
    path = []
    previous_node = None

    # loop over all the nodes in the user-specified line
    for loc in locs:
        # reshape the point
        test_pt = loc.reshape(-1, len(loc))

        # get the index and distance of the nearest node
        current_node, dist = nearest_node(verts, test_pt)

        # if this is not the first node
        if previous_node is not None:
            # find the shortest path between the previous node and this node
            djpath = nx.dijkstra_path(g, source=previous_node, target=current_node, weight='dist')
            # add the path the list of nodes (excluding the previous node)
            [path.append(x) for x in djpath[1:]]

        # if this is the first node
        else:
            # add it to the list
            path.append(current_node)

        # set the current node as the previous node
        previous_node = current_node

    return path


def get_string_by_card(all_strings, card=None, string_number=None):
    """Reduce the boundary_strings dataframe to a single card type or
    a single string number

    Args:
        all_strings: dataframe
            AdhSimulation.BoundaryCondition.boundary_strings dataframe
        card: string
            String card to extract from the dataframe (e.g. 'EGS', 'NDS', 'MTS')
        string_number: integer
            Number of the string to extract

    Returns:
        subset - dataframe
            subset dataframe containing only the requested cards and/or string

    """
    if card:
        card_condition = all_strings['CARD'] == card
    else:
        card_condition = all_strings['CARD'].notnull()

    if string_number:
        string_number_condition = all_strings['ID_1'] == string_number
    else:
        string_number_condition = all_strings['ID_1'].notnull()

    subset = all_strings.loc[card_condition & string_number_condition]

    return subset


def string_to_list(string_df, column_1='ID', column_2='ID_0'):
    """Convert AdH Boundary String dataframe into a list of properly ordered nodes
    AdH Boundary Strings are 2 column format, each row consists of two nodes that make
    up one edge. The second node is repeated as the first node on the next row

    Args:
        string_df: dataframe
            Boundary string dataframe (already reduced to a single string)
        column_1: string
            Dataframe column label of the first vertex of the edge
        column_2: string
            Dataframe column label for the second vertex of the edge

    Returns:
        Boundary string nodes as a list

    """
    # get the string series as a list
    lst = string_df[column_1].to_list()
    # add the last node in the string
    lst.append(string_df[column_2].iloc[-1])
    # return the list
    return lst


def get_boundary_string_path(mesh, node_list):
    """Convert a node list into path object

    Args:
        mesh: AdhMesh object
        node_list: List of nodes numbers with which to generate the path
            Nodes should be one-based since AdhMesh.verts index is one-based.

    Returns:
        Geoviews path object
    """

    path = gv.Path(mesh.verts.loc[node_list], crs=mesh.projection.get_crs())
    return path

