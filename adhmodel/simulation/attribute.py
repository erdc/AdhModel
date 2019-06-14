import networkx as nx
import math
from earthsim.annotators import GeoAnnotator
import param
from holoviews.streams import Selection1D, Stream
from holoviews import DynamicMap, Table
import pandas as pd

def attribute(fem_obj, start, stop, series=1):
    """
    Method to take a mesh and a start and stop location and map the edge path
    from one to the other.

    :param fem_obj: - obj
        roamsUTIL fem object
    :param start: - tuple
        Starting location of the edge (x1 float, y1 float)
    :param stop:  - tuple
        Ending location of the edge (x2 float, y2 float)
    :param series: - int
        Series number for this edge (default=1)
    :return:
    """
    # instantiate nx graph
    g = nx.Graph()

    # FEM to NetworkX
    # for each element in the mesh
    for element in fem_obj.elements:
        # add each element edge to the list of graph edges
        # (switch from one-based to zero based for indexing)
        g.add_edge(element[1] - 1, element[2] - 1, weight=0.0)
        g.add_edge(element[2] - 1, element[3] - 1, weight=0.0)
        g.add_edge(element[3] - 1, element[1] - 1, weight=0.0)

    # for each node in the mesh
    for node in fem_obj.nodes:
        # add to the graph nodes
        g.add_node(node[0] - 1, myIdx=node[0] - 1, xPos=node[1], yPos=node[2], zPos=node[3])

    # for every edge
    for u, v, d in g.edges(data=True):
        # calculate the length of the edge and store as 'weight'
        d['weight'] += math.sqrt((float(fem_obj.nodes[u][1]) - float(fem_obj.nodes[v][1])) ** 2 +
                                 (float(fem_obj.nodes[u][2]) - float(fem_obj.nodes[v][2])) ** 2)
    # loop over the nodes in the graph
    for node in g.nodes():
        # find the start node based on location  todo: add tolerance/rethink this
        if (float(g.nodes[node]['xPos']) == start[0]) & (float(g.nodes[node]['yPos']) == start[1]):
            startNodeIdx = int(g.nodes[node]['myIdx'])
        # find the stop node based on location  todo: add tolerance/rethink this
        elif (float(g.nodes[node]['xPos']) == stop[0]) & (float(g.nodes[node]['yPos']) == stop[1]):
            stopNodeIdx = int(g.nodes[node]['myIdx'])

    # find the shortest path between the start and stop nodes
    edge_path = nx.shortest_path(g, startNodeIdx, stopNodeIdx, d['weight'])
    # switch back to one-based indexing
    edge_path = [x + 1 for x in edge_path]
    # sort the answer
    path_sort = sorted(edge_path)

    print("BC")
    for edge in path_sort:
        print("NB ", edge, " ", series)


class AttributeAnnotator(GeoAnnotator):
    """
    Allows adding a group to the currently selected points using
    a dropdown menu and add button. The current annotations are
    reflected in a table.

    Works by using a shared datasource on the Points and Table.
    """

    add = param.Action(default=lambda self: self.add_group(), precedence=1,
                       doc="""Button triggering add action.""")

    group = param.ObjectSelector()

    column = param.String(default='Group', constant=True)

    table_height = param.Integer(default=150, doc="Height of the table",
                                 precedence=-1)

    def __init__(self, groups, **params):
        super(AttributeAnnotator, self).__init__(**params)
        self.param.group.objects = groups
        self.param.group.default = groups[0]
        # self.point_sel_stream = Selection1D(source=self.points)
        self.poly_sel_stream = Selection1D(source=self.polys)
        self._group_data = {g: [] for g in groups}
        self.table_stream = Stream.define('TableUpdate')(transient=True)

    def add_group(self, **kwargs):
        new_index = self.point_sel_stream.index
        for idx in new_index:
            if idx not in self._group_data[self.group]:
                self._group_data[self.group].append(idx)
            for g, inds in self._group_data.items():
                if g != self.group:
                    self._group_data[g] = [idx for idx in inds if idx not in new_index]
        self.table_stream.trigger([self.table_stream])

    def group_table(self):
        plot = dict(width=self.width, height=self.table_height)
        data = [(group, str(inds)) for group, inds in self._group_data.items()]
        return Table(data, self.column, 'index').sort().opts(plot=plot)

    def annotated_polys(self):
        element = self.poly_stream.element
        groups = []
        for g, idx in self._group_data.items():
            df = element.iloc[idx].dframe()
            df[self.column] = g
            groups.append(df)
        data = pd.concat(groups).sort_values(self.column) if groups else []
        return element.clone(data, vdims=self.column).opts(plot={'color_index': self.column},
                                                           style={'cmap': 'Category20'})

    def view(self):
        table = DynamicMap(self.group_table, streams=[self.table_stream])
        annotated = DynamicMap(self.annotated_polys, streams=[self.table_stream])
        return (self.tiles * self.polys * self.points * annotated + table).cols(1)



# path = 'C:/PROJECTS/ERS'
# filename = 'angle_db'
# start = (0, 0)
# stop = (21.2132034356, -21.2132034356)
#
# adh_obj = model.adhModel()
# fem_obj = adh_obj.readMeshFile(path, filename)
#
# attribute(fem_obj, start, stop)



