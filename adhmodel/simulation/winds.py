import datetime
from datetime import timedelta
import os
import pandas as pd
from itertools import chain
from scipy.spatial import distance
import param
import cartopy.crs as ccrs
import numpy as np

# from bokeh.layouts import gridplot
# from bokeh.plotting import figure, show, output_file
# from bokeh.io import export_png
import holoviews as hv
import geoviews as gv


class OWI(param.Parameterized):
    grid = None
    field1 = pd.DataFrame()
    field2 = pd.DataFrame()
    delta = None
    number_fields = param.Integer(2, bounds=(1, 2), precedence=2)
    ilat = param.Integer(0, precedence=2)
    ilon = param.Integer(0, precedence=2)
    dx = param.Number(0.25, precedence=2)
    dy = param.Number(0.25, precedence=2)
    swlat = param.Number(5.8, precedence=2)
    swlon = param.Number(5.8, precedence=2)
    time_range = param.DateRange(default=None, doc="""
                    Start and end datetime""")
    date_start = param.Date(default=datetime.datetime.now(), precedence=-1)
    date_stop = param.Date(default=datetime.datetime.now() + timedelta(hours=5), precedence=-1)

    def get_file_metadata(self, path, filename):
        """
        Method to get extract only time between snaps from an OWI file using
        the first two snaps.

        Parameters
        ----------
        path - str
            directory where the file is located
        filename - str
            filename of the owi file

        Returns
        -------
        delta - datetime.timedelta
            time between snaps
        """
        # construct full path
        full_path = os.path.join(path, filename)
        # open the file
        fp = open(full_path)
        # get the first line
        line1 = fp.readline().split()
        # date_start = datetime.datetime.strptime(line1[-2], '%Y%m%d%H')
        # date_stop = datetime.datetime.strptime(line1[-1], '%Y%m%d%H')
        # get the second line
        line2 = fp.readline().rstrip('\n')
        # get the first snap time
        start_time = datetime.datetime.strptime(line2[68:80], '%Y%m%d%H')

        for line in fp:
            # if this is a header line
            if line[0:2] == 'iL':
                line.rstrip('\n')
                # get the second snap time
                snap_time = datetime.datetime.strptime(line.rstrip('\n')[68:80], '%Y%m%d%H')
                # calculate the time between snaps
                delta = snap_time - start_time

                fp.close()

                return delta

        fp.close()

        print('An accurate delta could not be calculated because only one snap '
              'was found in the file. Default of 0.0 will be used')

        delta = datetime.timedelta(seconds=0.0)
        return delta

    def read_owi_wind(self, path, filename, fields):
        """ Method to read owi wind field from file

        Parameters
        ----------
        path - str
            Directory where the file is located
        filename - str
            Filename to be read (with extension)
        delta - datetime.timedelta
            time between time snaps
        fields - int
            Number of fields in this file (1 for single data, 2 for x and y data)

        """
        if self.delta is None:
            self.delta = self.get_file_metadata(path, filename)

        # construct full path
        full_path = os.path.join(path, filename)

        # # convert the provided delta to a timedelta interval
        # interval = datetime.timedelta(minutes=delta)

        # set the nodes per line (fixed OWI format)
        nodes_per_line = 8

        # instantiate variables
        snaps = []
        data = []
        next_header = 0
        header_index = 0
        field1_df = pd.DataFrame()
        field2_df = pd.DataFrame()

        # open the file
        fp = open(full_path)
        # loop over all the lines in the file (with an index)
        for i, line in enumerate(fp):
            # first line has file header information
            if i == 0:
                # parse the line
                spar = line.split()
                # record the start time
                time_start = datetime.datetime.strptime(spar[3], '%Y%m%d%H')
                # record the stop time
                time_stop = datetime.datetime.strptime(spar[4], '%Y%m%d%H')
                # calculate the number of snaps in the file
                (number_snaps, remainder) = divmod((time_stop - time_start).total_seconds(), self.delta.total_seconds())
                # check for even division
                if remainder != 0.0:
                    raise UserWarning('Warning: The time delta did not evenly match with the start and '
                                      'stop time in the file. The results may contain bad data.')

                # calculate the datetime of each snap
                snaps = [time_start + self.delta * x for x in range(0, int(number_snaps) + 1, 1)]

            # the first data header line
            elif i == 1:
                # parse the line
                # get number of latitude nodes
                ilat = int(line[5:9])
                # get number of longitude nodes
                ilon = int(line[15:19])
                # get dx
                dx = float(line[22:28])
                # get dy
                dy = float(line[31:37])
                # get the SW corner of the grid latitude
                swlat = float(line[43:51])
                # get the sW corner of the grid longitude
                swlon = float(line[57:65])
                # get the snap time as a string
                snap_time_str = line[68:80]

                # calculate the number of nodes
                number_nodes = ilat * ilon
                # calculate the number of lines per snap (per field)
                lines_per_snap = int(np.ceil(number_nodes / nodes_per_line))
                # calculate the indices for the header lines
                header_lines = [x for x in range(1, 1 + (len(snaps) * (1 + lines_per_snap * fields)),
                                                 (1 + lines_per_snap * fields))]

                # increment the header line index
                header_index += 1
                # get the next header line
                next_header = header_lines[header_index]

            # stop at the next header line
            elif i == next_header:
                # flatten the data list
                column_data = list(chain.from_iterable(data))
                if fields == 1:
                    # store the previous data for field1
                    field1_df[snap_time_str] = column_data[0:len(column_data)]

                # if there are two fields in the file
                else:
                    # store the previous data for field1
                    field1_df[snap_time_str] = column_data[0:int(len(column_data) / 2)]
                    # store the field2 data
                    field2_df[snap_time_str] = column_data[int(len(column_data) / 2):len(column_data)]


                # record the snap time for the next interval
                snap_time_str = line[68:80]

                # clear the data list
                data = []

                # increment the header flag
                header_index += 1
                # if we've just moved beyond the last header in the list, stop
                if header_index == len(header_lines):
                    # break
                    next_header = -1
                else:
                    # get the next header line
                    next_header = header_lines[header_index]

            # if this is not a header line, it is data, append to the data list
            else:
                data.append([float(x) for x in line.split()])

        # Store the data from the last timestep
        # flatten the data list
        column_data = list(chain.from_iterable(data))
        if fields == 1:
            # store the previous data for field1
            field1_df[snap_time_str] = column_data[0:len(column_data)]

        # if there are two fields in the file
        else:
            # store the previous data for field1
            field1_df[snap_time_str] = column_data[0:int(len(column_data) / 2)]
            # store the field2 data
            field2_df[snap_time_str] = column_data[int(len(column_data) / 2):len(column_data)]

        # close the file
        fp.close()

        # store data in the object
        # self.delta = delta
        self.number_fields = fields
        self.ilat = ilat
        self.ilon = ilon
        self.dx = dx
        self.dy = dy
        self.swlat = swlat
        self.swlon = swlon
        # self.start_time = time_start
        # self.stop_time = time_stop
        self.time_range = (time_start, time_stop) # todo

        # if there are two fields, return both
        if fields == 2:
            self.field1 = field1_df
            self.field2 = field2_df
        # otherwise return field1
        else:
            self.field1 = field1_df

    def set_grid_parameters(self, swLat, swLon, neLat, neLon, dx, dy):
        """Method to set the actual grid values based on a bounding box.
        dx and dy default to 0.25. If the number of cells cannot be
        evenly divided, an extra grid cell is added. Values are set into self.
        Latitude = y
        Longitude = x
        """

        # calculate the estimated grid size
        delta_lat = abs(swLat - neLat)
        delta_lon = abs(swLon - neLon)

        # calculate the number of cells required
        ilat, rem = divmod(delta_lat, dy)
        # if the delta cannot be evenly divided, add a cell
        if rem != 0:
            ilat += 1
        # if the ilat is even, add a cell
        if ilat % 2 == 0:
            ilat += 1
        # calculate the number of cells required
        ilon, rem = divmod(delta_lon, dx)
        # if the delta cannot be evenly divided, add a cell
        if rem != 0:
            ilon += 1
        # if the ilon is even, add a cell
        if ilon % 2 == 0:
            ilon += 1

        # set values into self
        self.ilat = int(ilat)
        self.ilon = int(ilon)
        self.dx = dx
        self.dy = dy
        self.swlat = swLat
        self.swlon = swLon

    def make_owi_grid(self):
        """
        Method to construct the list of grid nodes (latitude and longitude) from
        the construction parameters. The indexes of these grid nodes correspond to
        the indexing of the data in the OWI file.

        Parameters
        ----------
        swlon - float
            Longitude of the SouthWestern Corner of the grid
        swlat - float
            Latitude of the SouthWestern Corner of the grid
        dx - float
            grid spacing in the x direction (Longitude)
        dy - float
            grid spacing in the y direction (Latitude)
        ilon - int
            number of grid cells in the x direction (Longitude)
        ilat - int
            number of grid cells in the y direction (Latitude)

        Returns
        -------
        self.grid - df
            DataFrame containing the Latitude and Longitude lists for all the grid
            points. The indexing corresponds to the indexing of the data arrays from
            the file.
        """
        # ensure that all the required values have been set
        if self.ilat is None or self.ilon is None or self.dx is None or self.dy is None or \
                        self.swlat is None or self.swlon is None:
            raise IOError('Parmeters defining the grid must be set into the OWI object (ilat, '
                          'ilon, dx, dy, swlat, swlon)')

        # create array of x values
        xs = np.arange(self.swlon, self.swlon + self.ilon * self.dx, self.dx)

        # create array of y values
        ys = np.arange(self.swlat, self.swlat + self.ilat * self.dy, self.dy)

        # construct grid
        data = np.stack([[x, y] for y in ys for x in xs], axis=1)

        # construct the dataframe for the grid
        self.grid = pd.DataFrame({'Longitude': data[0], 'Latitude': data[1]})

        # ensure all grid values are valid
        self.check_grid()

    def check_grid(self):
        """Method to check that the grid does not cross a pole or an
        antimeridian with invalid values. Corrects as needed and sets
        values back into self.grid."""

        if self.grid.Longitude.min() < -180.0:
            bool = self.grid.Longitude < -180.0
            self.grid.Longitude[bool] = self.grid.Longitude[bool] + 180.0

        if self.grid.Longitude.max() > 180.0:
            bool = self.grid.Longitude > 180.0
            self.grid.Longitude[bool] = self.grid.Longitude[bool] - 180

        if self.grid.Latitude.min() < -90.0:
            bool = self.grid.Latitude < -90.0
            self.grid.Latitude[bool] = self.grid.Latitude[bool] + 90.0

        if self.grid.Latitude.min() > 90.0:
            bool = self.grid.Latitude > 90.0
            self.grid.Latitude[bool] = self.grid.Latitude[bool] - 90.0

    def create_wind_field(self, wind_df, directory, filename):
        """Method to generate wind field. Derived from Holland 1980"""
        #ensure data is present
        if wind_df.empty:
            raise IOError('No data provided. Unable to create wind field.')

        # Create the grid for this wind field (from the dx, dy, ilat, ilon, swlat, swlon)
        if self.grid is None:
            self.make_owi_grid()

        for index in wind_df.index:
            # date = wind_df['Datetime'][index]
            # lat = wind_df['Latitude'][index]
            # lon = wind_df['Longitude'][index]
            # sf_kt = wind_df['sf_kt'][index]
            # Rw_nm = wind_df['Rw_nm'][index]
            # pc_mb = wind_df['pc_mb'][index]
            # deg = wind_df['deg'][index]
            # v_kt = wind_df['v_kt'][index]

            date = wind_df['Datetime'][index]
            lat = wind_df['Latitude'][index]
            lon = wind_df['Longitude'][index]
            sf_kt = wind_df['Max_Wind_Speed'][index]
            Rw_nm = wind_df['Radius'][index]
            pc_mb = wind_df['Pressure'][index]
            deg = wind_df['Direction_Degrees'][index]
            v_kt = wind_df['Direction_Speed'][index]

            # convert grid to np array format # todo there is probably a better way to do this
            gridarray = np.array(list(zip(self.grid.Longitude, self.grid.Latitude)))

            # set the hurricane center point
            pt = np.array([lon, lat])
            # force the second axis for matrix multiplication
            pt = pt[np.newaxis, ...]
            # calculate the distances from the hurricane point
            degrees = distance.cdist(gridarray, pt)
            degrees = degrees.reshape(len(degrees))
            # set the avg radius of the earth
            earth_radius = 6371 * 1000  # meters
            # convert degrees to radius
            r_array = (np.pi * earth_radius) * degrees / 180
            delta_lon = (np.pi * earth_radius) * (self.grid.Longitude.values - lon) / 180 # todo check this
            delta_lat = (np.pi * earth_radius) * (self.grid.Latitude.values - lat) / 180
            # ignore div by zero errors
            np.seterr(divide='ignore')
            # calculate thea
            theta_array = np.arctan2(delta_lat, delta_lon)

            #######################################
            # set standard assumptions
            # boundary layer adjustment factor
            beta = 0.9
            # density of air
            den_air = 1.15  # kgm^-3
            # density of water
            den_water = 997.  # kgm^-3
            # ambient atmospheric pressure https://www.britannica.com/science/atmospheric-pressure
            pn = 101325
            # Euler's number
            e = 2.718281
            # gravity
            grav = 9.80665  # ms^2
            # rotation rate of the earth
            omega = 7.2921e-5  # rad/s (2*pi/day)
            # sampling adjustment (converts 1 min winds to 10 min winds)
            ct = 0.88
            ########################################################
            # apply conversions
            # convert knots to mps
            sf = sf_kt * 0.5144444
            v = v_kt * 0.5144444
            # convert nm to meters
            Rw = Rw_nm * 1852
            # convert millibars to pascal
            pc = pc_mb * 100
            # convert translational speed to
            # velocity east
            vte = v * np.sin(deg * np.pi / 180)
            # velocity west
            vtn = v * np.cos(deg * np.pi / 180)
            ######################################################

            # calculate maximum storm wind speed at 10m
            sm = sf - np.sqrt(np.power(vte, 2) + np.power(vtn, 2))
            # calculate maximum velocity at the top of the atmospheric boundary layer
            Vm = sm / beta
            # calculate Holland's B parameter
            holland_b = den_air * e * np.power(Vm, 2) / (pn - pc)
            # ensure Holland's B parameter is within reasonable limits 1 < B < 2.5
            if holland_b < 1:
                # set to lower limit of 1
                holland_b = 1
            elif holland_b > 2.5:
                # set to upper limit of 2.5
                holland_b = 2.5

            # coriolis parameter for the storm's current location
            f = 2 * omega * np.sin(lat * np.pi / 180)  # todo replace with array # typical values about 10^-4 rads/s

            # for each node from the center of storm and it's radial angle

            # calculate pressure at each node
            # pr_array = (pc + (pn - pc) * np.exp(-np.power(Rw / r_array, holland_b))) / (den_water * grav)

            # calculate the raw gradient wid speed at each node
            # Vg = np.sqrt((np.power(Rw / r, holland_b) - np.exp(1 - np.power(Rw / r, holland_b)) * np.power(Vm, 2) +
            #               (np.power(float(r), 2) * np.power(f, 2) / 4))
            #              ) g- (r * f / 2) # orig - suspected typo
            Vg = np.sqrt(
                (np.power(Rw / r_array, holland_b) * np.exp(1 - np.power(Rw / r_array, holland_b)) * np.power(Vm, 2) +
                 (np.power(r_array, 2) * np.power(f, 2) / 4))) - (r_array * f / 2)  # derived from Holland 1980

            # translational adjustments to final wind speed in north and east directions
            vtan = np.abs(Vg / Vm) * vtn
            vtae = np.abs(Vg / Vm) * vte

            # for each node location, i, separate the velocity into components and apply beta
            Vei = -Vg * beta * np.sin(theta_array)  # todo: does this need to be adjusted for southern hemi?
            Vni = Vg * beta * np.cos(theta_array)

            # multiply by sampling time adjustment to convert 1min winds to 10min winds and add translation velocity
            Vfei = ct * Vei + vtae
            Vfni = ct * Vni + vtan

            date_string = datetime.datetime.strftime(date, '%Y%m%d%H%M')
            self.field1[date_string] = Vfei
            self.field2[date_string] = Vfni

        # set start and stop dates
        self.date_start = wind_df['Datetime'].min()
        self.date_stop = wind_df['Datetime'].max()

        self.write_owi_file(directory, filename)

    def write_owi_file(self, path, filename):
        """ Method to write owi formatted wind file from values in the owi object.

        Parameters
        ----------
        path - str
          Output directory
        filename - str
          Output filename (with extension)

        Returns
        -------


        """
        # construct full path
        full_path = os.path.join(path, filename)

        # ensure that all the required values have been set
        if self.ilat is None or self.ilon is None or self.dx is None or self.dy is None or \
                        self.swlat is None or self.swlon is None:
            raise IOError('Parmeters defining the grid must be set into the OWI object (ilat, '
                          'ilon, dx, dy, swlat, swlon)')
        if self.grid is None:
            raise IOError('Grid must be constructed and set in the OWI object.')
        if self.field1 is None:
            raise IOError('At least one field of data must be set in the OWI object')
        if self.number_fields == 2 and self.field2 is None:
            raise IOError('Number of fields specified as 2, but no field2 data was found in the OWI object')
        # if self.start_time is None and self.stop_time is None:
        #     raise IOError('Start time and stop time must be specified as %Y%m%d%H')
        # set the number of columns in the file (fixed OWI format)
        columns = 8

        # start_time = self.time_range[0]
        # stop_time = self.time_range[1]
        start_time = self.date_start  # todo
        stop_time = self.date_stop  # todo

        # open the file
        fp = open(full_path, 'w+', newline='')

        # construct file header
        file_header = 'Oceanweather WIN/PRE Format                            {0}     {1}\n'.\
            format(datetime.datetime.strftime(start_time, '%Y%m%d%H'),
                   datetime.datetime.strftime(stop_time, '%Y%m%d%H'))
        # write the file header
        fp.write(file_header)

        # loop over the timesteps in the data
        for time_string in self.field1.__iter__():
            # construct the timestep header
            timestep_header = 'iLat={:>4}iLong={:>4}DX={:01.4f}DY={:01.4f}SWLat={:08.5f}SWLon={:02.4f}DT={}\n'.format(
                self.ilat, self.ilon, self.dx, self.dy, self.swlat, self.swlon, time_string)
            # write the timestep header
            fp.write(timestep_header)

            # loop over the data fields
            for data in [self.field1[time_string].values, self.field2[time_string].values] if self.number_fields == 2 \
                    else [self.field1[time_string].values]:

                # determine the number of rows required per field and the remaining extra data columns
                (num_rows, rem_columns) = divmod(len(data), columns)

                # reformat into fixed width strings and construct the columns for the file
                fullrows = [['{:>10.4f}'.format(data[i + j]) for j in range(columns)] for i in
                            range(0, len(data) - rem_columns, 8)]
                # print each full row
                for row in fullrows:
                    row = ''.join(row)
                    fp.write(row + '\n')
                # print the last row of data (doesn't contain all the columns)
                lastrow = ['{:>10.4f}'.format(data[num_rows * 8 + j]) for j in range(rem_columns)]
                fp.write(''.join(lastrow) + '\n')

    def write_grid_xy(self, path, filename):
        """ Method to write owi grid file as XY from values in the owi object.

        Parameters
        ----------
        path - str
          Output directory
        filename - str
          Output filename (with extension)

        Returns
        -------


        """
        # construct full path
        full_path = os.path.join(path, filename)

        if self.grid is None:
            raise IOError('Grid must be constructed and set in the OWI object.')

        self.grid.to_csv(full_path, sep='\t', index=False)

    def view(self):

        def time_field(time):
            vx = self.field1[time].values
            vy = self.field2[time].values
            xs = self.grid.Longitude.values
            ys = self.grid.Latitude.values
            with np.errstate(divide='ignore', invalid='ignore'):
                angle = np.arctan2(vy, vx)
            mag = np.sqrt(vx ** 2 + vy ** 2)
            return gv.VectorField((xs, ys, angle, mag), vdims=['Angle', 'Magnitude'],
                                  crs=ccrs.PlateCarree()).redim.range(Magnitude=(0, 50))

        # tiles = gv.WMTS('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{Z}/{Y}/{X}')
        vectors = hv.DynamicMap(time_field, kdims='Time').redim.values(Time=sorted(self.field1.keys()))

        return vectors.opts(size_index='Magnitude', color_index='Magnitude', pivot='tail',
                            width=700, height=500, colorbar=True, scale=1, cmap='OrRd')
