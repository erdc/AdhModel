import unittest
import os
import filecmp
import pandas as pd
from adhmodel.adh_model import AdhModel
from adhmodel.simulation.simulation import AdhSimulation, BoundaryConditions
from adhmodel.simulation.io import read_bc_file, write_bc_file, write_hot_start_file
import holoviews.plotting.bokeh
import geoviews.plotting.bokeh


class TestIo(unittest.TestCase):

    def test_read_mesh_file(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_bt.3dm')
        out_file = os.path.join(path, 'out_files', 'out_angle_bt.3dm')
        base_file = os.path.join(path, 'base_angle_bt2.3dm')
        model = AdhModel()
        model.read_mesh(fname, fmt='3dm')
        model.write_mesh(out_file, fmt='3dm')
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_bad_format_nds(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'bad_mesh_nds.3dm')
        model = AdhModel()
        with self.assertRaises(UnboundLocalError) as context:
            model.read_mesh(fname, fmt='3dm')
        # self.assertTrue('Error reading node (ND) record from file on line 2.' in str(context.exception))
        self.assertTrue(
            "local variable 'extra_cols' referenced before assignment" in str(context.exception))

    def test_bad_format_e3t(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'bad_mesh_e3t.3dm')
        model = AdhModel()
        with self.assertRaises(ValueError) as context:
            model.read_mesh(fname, fmt='3dm')
        # self.assertTrue('Error reading element (E3T) record from file on line 3.' in str(context.exception))
        self.assertTrue("list.remove(x): x not in list" in str(context.exception))

    def test_read_hot_start(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_bt.hot')
        out_file = os.path.join(path, 'out_files', 'out_angle_bt.hot')
        base_file = os.path.join(path, 'base_angle_bt.hot')
        model = AdhModel()
        model.read_hotstart(fname, fmt='ascii')
        hotstart = model.simulation.hotstart
        self.assertEqual(1, len(hotstart))
        self.assertEqual(250, hotstart.ioh.size)
        model.simulation.write_hotstart(out_file)
        self.assertTrue(filecmp.cmp(base_file, out_file, shallow=False))

    def test_read_hot_start_2(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'ka_mult_sed.hot')
        out_file = os.path.join(path, 'out_files', 'out_ka_mult_sed.hot')
        base_file = os.path.join(path, 'base_ka_mult_sed.hot')
        model = AdhModel()
        model.read_hotstart(fname, fmt='ascii')
        hotstart = model.simulation.hotstart
        self.assertEqual(2, len(hotstart))
        self.assertEqual(10891, hotstart.ioh.size)
        self.assertEqual((1, 10891, 3), hotstart.iov.shape)
        model.simulation.write_hotstart(out_file)
        self.assertTrue(filecmp.cmp(base_file, out_file, shallow=False))

    def test_read_bc_red_river(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'RedRiver.bc')
        base_file = os.path.join(path, 'base_RedRiver.bc')
        out_file = os.path.join(path, 'out_files', 'out_RedRiver.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_red_river2(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'RR_Fine_12052017_1038.bc')
        base_file = os.path.join(path, 'base_RR_Fine_12052017_1038.bc')
        out_file = os.path.join(path, 'out_files', 'out_RR_Fine_12052017_1038.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_east_ll_coarse(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'east_LL_coarse.bc')
        base_file = os.path.join(path, 'base_east_LL_coarse.bc')
        out_file = os.path.join(path, 'out_files', 'out_east_LL_coarse.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_anc_UTM(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'anc_UTM.bc')
        base_file = os.path.join(path, 'base_anc_UTM.bc')
        out_file = os.path.join(path, 'out_files', 'out_anc_UTM.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_hunter(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'hunter.bc')
        base_file = os.path.join(path, 'base_hunter.bc')
        out_file = os.path.join(path, 'out_files', 'out_hunter.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_2d_coriolis_x(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, '2D-coriolis_X.bc')
        base_file = os.path.join(path, 'base_2D-coriolis_X.bc')
        out_file = os.path.join(path, 'out_files', 'out_2D-coriolis_X.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_angle_dam(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_dam.bc')
        base_file = os.path.join(path, 'base_angle_dam.bc')
        out_file = os.path.join(path, 'out_files', 'out_angle_dam.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_angle_nb(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_nb.bc')
        base_file = os.path.join(path, 'base_angle_nb.bc')
        out_file = os.path.join(path, 'out_files', 'out_angle_nb.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_angle_source(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_source.bc')
        base_file = os.path.join(path, 'base_angle_source.bc')
        out_file = os.path.join(path, 'out_files', 'out_angle_source.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_angle_wave(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_wave.bc')
        base_file = os.path.join(path, 'base_angle_wave.bc')
        out_file = os.path.join(path, 'out_files', 'out_angle_wave.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_angle_wind(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'angle_wind.bc')
        base_file = os.path.join(path, 'base_angle_wind.bc')
        out_file = os.path.join(path, 'out_files', 'out_angle_wind.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_wd(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'wd.bc')
        base_file = os.path.join(path, 'base_wd.bc')
        out_file = os.path.join(path, 'out_files', 'out_wd.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

    def test_read_bc_file(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'ka_mult_sed.bc')
        base_file = os.path.join(path, 'base_ka_mult_sed.bc')
        out_file = os.path.join(path, 'out_files', 'out_ka_mult_sed.bc')
        out_file2 = os.path.join(path, 'out_files', 'out2_ka_mult_sed.bc')
        model = AdhModel()
        model.read_bc(fname, fmt='bc')
        model.write_bc(out_file, fmt='bc', validate=False)
        self.assertTrue(filecmp.cmp(base_file, out_file))

        bc_class2 = BoundaryConditions()
        read_bc_file(out_file, bc_class2)
        write_bc_file(out_file2, bc_class2)
        self.assertTrue(filecmp.cmp(out_file, out_file2))

    def test_read_hot_start_bad_file(self):
        pass
        # todo commented out because the hot reader from nmutils does not read line by line
        # todo and would therefore require all of these tests to be refactored.

        # path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        # fname = os.path.join(path, 'out_files', 'bad.hot')

        # file_str = 'bad first line\n\n'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue('First line in file must be "DATASET".' in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nBEGSCL\n'
        # msg = 'Invalid format on line 3. Found "BEGSCL" but currently reading data set.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nBEGVEC\n'
        # msg = 'Invalid format on line 3. Found "BEGVEC" but currently reading data set.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nND\nBEGSCL\n'
        # msg = 'Invalid format on line 2. Found "ND" but no "BEGSCL" or "BEGVEC" previously read.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nND\n'
        # msg = 'Invalid format on line 3. Unable to read number of nodes (ND).'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nNC\nBEGSCL\n'
        # msg = 'Invalid format on line 2. Found "NC" but no "BEGSCL" or "BEGVEC" previously read.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nNC\n'
        # msg = 'Invalid format on line 3. Unable to read number of cells (NC).'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nNAME\nBEGSCL\n'
        # msg = 'Invalid format on line 2. Found "NAME" but no "BEGSCL" or "BEGVEC" previously read.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nNAME\n'
        # msg = 'Invalid format on line 3. Unable to read name.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nNAME ""\n'
        # msg = 'Invalid format on line 3. Unable to read name.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nTS\nBEGSCL\n'
        # msg = 'Invalid format on line 2. Found "TS" but no "BEGSCL" or "BEGVEC" previously read.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nTS\n'
        # msg = 'Invalid format on line 3. Unable to read time index and time (TS).'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nTS 0\n'
        # msg = 'Invalid format on line 3. Unable to read time index and time (TS).'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nENDDS\nBEGSCL\n'
        # msg = 'Invalid format on line 2. Found "ENDDS" but no "BEGSCL" previously read.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\nND 3\nNC 1\n1\n2\nENDDS\n'
        # msg = 'Invalid format on line 7. Found "ENDDS" but number of data set values read does not equal number ' \
        #       'of nodes (ND).'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = 'DATASET\nBEGSCL\n1\n'
        # msg = 'Invalid format on line 3. Number of nodes not read (ND).'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = '''DATASET
        #             OBJTYPE "mesh2d"
        #             BEGSCL
        #             ND 8
        #             NC 8
        #             NAME "sediment transport"
        #             RT_JULIAN 2453867.068720
        #             TIMEUNITS seconds
        #             TS 1 1.00000000e+00
        #             0
        #             0
        #             0
        #             1
        #             1
        #             1
        #             1
        #             __BAD__LINE__
        #             0.00000000e+00
        #             0.00000000e+00
        #             0.00000000e+00
        #             3.24000000e+00
        #             4.39000000e+00
        #             2.96000000e+00
        #             7.48000000e+00
        #             0.00000000e+00
        #             ENDDS
        #             '''
        # msg = 'Invalid format on line 17. Unable to read data set activity value.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

        # file_str = '''DATASET
        #             OBJTYPE "mesh2d"
        #             BEGSCL
        #             ND 8
        #             NC 8
        #             NAME "sediment transport"
        #             RT_JULIAN 2453867.068720
        #             TIMEUNITS seconds
        #             TS 1 1.00000000e+00
        #             0
        #             0
        #             0
        #             1
        #             1
        #             1
        #             1
        #             0
        #             0.00000000e+00
        #             0.00000000e+00
        #             0.00000000e+00
        #             3.24000000e+00
        #             4.39000000e+00
        #             2.96000000e+00
        #             7.48000000e+00
        #             __BAD__LINE__
        #             ENDDS
        #             '''
        # msg = 'Invalid format on line 25. Unable to read data set value.'
        # with open(fname, "w") as hot_file:
        #     hot_file.write(file_str)
        # with self.assertRaises(IOError) as context:
        #     model = AdhModel()
        #     model.read_hotstart(fname, fmt='3dm')
        # self.assertTrue(msg in str(context.exception))

    def test_read_bc_bad_file(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'out_files', 'bad.bc')
        bc_class = BoundaryConditions()

        file_str = 'OP TRN\n'
        msg = 'Error reading line 1 of file: bad.bc.\nLine: OP TRN'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

    def test_read_bc_misc_cards(self):
        file_str = '''OP SW3
                    OP BT
                    OP BTS
                    OP TEM 0.65
                    OP TPG 0.45
                    OP NF2
                    OP WND
                    OP WAV
                    OP DAM
                    OP DIF
                    OP TRN 4

                    IP NTL 0.35
                    IP ITL 0.44

                    CN SAL 1 19.0
                    CN TMP 2 25.0 1
                    CN VOR 3 1.0 2.0 3.0
                    CN CON 4 .03

                    MP WND STR 1 0
                    MP WND ATT 1 0.5
                    MP EEV 1 1.23 4
                    MP COR 1 45.0

                    FR MNG 1 .03
                    FR MNC 2 .035
                    FR ERH 1 1.23
                    FR SAV 1 2.34 3.45
                    FR URV 2 4.56 5.67 6.78
                    FR EDO 3 1 2 3 4 5
                    FR ICE 1 9.87 8.76 1
                    FR IRH 1 .03
                    FR BRH 1 3.21
                    FR SDK 1 4.32
                    FR BRD 2 9.99 3.33

                    FLX 1
                    FLX 2
                    PC ADP
                    PC ELM
                    PC MEO

                    DB RAD 3 4
                    NB SDR 5 5 4 3 2 1
                    BR JAI 1 0 44.4 33.3 22.2 11.1 9 8
                    BR SAS 1 0 99 88 77 66 2 3
                    BR MLM 999 888 777 666
                    BR FRO 1 0 11 22 33 44 55 4 5 
                    BR BRC 1 0 11 22 33 44 4 5 
                    BR VTG 1 0 44.4 33.3 22.2 1 55.5 6 7
                    BR FER 1 0 33.3 22.2 11.1 0 44.4 9 8 
                    BR USR 1 2
                    WER 1
                    WRS 1 3 4 5 6 55.5 66.6 77.7
                    FLP 1
                    FGT 1 2 5 4 7 6 1.1 2.2 3.3 4.4 5.5 6.6 77.7
                    SLUICE 1
                    SLS 1 5 6 9 8 55.5 2
                    OFF 5
                    OB OF 4


                    SERIES WIND 1 2 1.5 2.6 0 0
                    1 2 3
                    4 5 6

                    SERIES WAVE 2 2 1.5 2.6 0 0
                    11 12 13
                    14 15 16

                    TC T0 0 3
                    TC ATF 2.2 3
                    TC STD 3.3 4.4

                    PC LVL 1
                    PC LVL 0
                    SOUT RESID
                    SOUT MERROR
                    SOUT NLNODE
                    SOUT LNODE
                    FOUT WIND
                    FOUT WAVE
                    FOUT ADAPT GRID
                    FOUT ADAPT SW
                    FOUT ADAPT CON
                    FOUT SED
                    '''
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'out_files', 'misc_cards.bc')
        out_fname = os.path.join(path, 'out_files', 'out_misc_cards.bc')
        base_file = os.path.join(path, 'base_misc_cards.bc')
        bc_class = BoundaryConditions()
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        read_bc_file(fname, bc_class)
        write_bc_file(out_fname, bc_class)
        self.assertTrue(filecmp.cmp(base_file, out_fname))

    def test_read_bc_with_errors(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'out_files', 'misc_cards.bc')
        out_fname = os.path.join(path, 'out_files', 'out_2.bc')
        bc_class = BoundaryConditions()

        file_str = '''IP ITL .22\nIP NTL .35\nIP NTL bob\n'''
        msg = 'Error reading line 3 of file: misc_cards.bc.\nLine: IP NTL bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''CN CON bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: CN CON bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''MP WND STR bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: MP WND STR bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''XY1 bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: XY1 bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''NDS bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: NDS bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''FR MNG bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: FR MNG bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''PC LVL bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: PC LVL bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''NB DIS bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: NB DIS bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''WRS bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: WRS bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''TC T0 bob\n'''
        msg = 'Error reading line 1 of file: misc_cards.bc.\nLine: TC T0 bob'
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        with self.assertRaises(IOError) as context:
            read_bc_file(fname, bc_class)
        self.assertEqual(msg, str(context.exception))

        file_str = '''TC ATF 2.2 3\n'''
        with open(fname, "w") as bc_file:
            bc_file.write(file_str)
        read_bc_file(fname, bc_class)
        write_bc_file(out_fname, bc_class)
        with open(out_fname, "r") as bc_file:
            tmp_str = '\n'.join(bc_file.readlines())
            self.assertTrue('TC ATF 2.2 3' in tmp_str)

        recs = [['BR', 'JAI', 1, 1, 1.1, 2.2, 3.3, 4.4, 1, 1],
                ['BR', 'SAS', 1, 1, 1.1, 2.2, 3.3, 4.4, 1, 1], ]
        bc_class.breach_controls = pd.DataFrame.from_records(recs)
        msg = 'Error writing BR card.'
        with self.assertRaises(IOError) as context:
            write_bc_file(out_fname, bc_class)
        self.assertEqual(msg, str(context.exception))

    def test_validate_parameters(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files')
        fname = os.path.join(path, 'base_anc_UTM.bc')
        bc_class = BoundaryConditions()
        read_bc_file(fname, bc_class)
        msg = bc_class.validate_sw2()
        self.assertTrue(msg)

    def test_read_results(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files/SanDiego')
        project_name = 'SanDiego'
        model = AdhModel()
        model.read_results(path=path, project_name=project_name, fmt='ascii')
        results = model.simulation.results
        self.assertEqual(len(results), 3)
        self.assertIn('Depth', results)
        self.assertIn('Error', results)
        self.assertIn('Depth-Averaged Velocity', results)
        self.assertEqual(results.Depth.shape, (49, 9140))
        self.assertEqual(results['Depth-Averaged Velocity'].shape, (49, 9140, 3))

    def test_read_result(self):
        file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_files/SanDiego', 'SanDiego_dep.dat')
        model = AdhModel()
        model.read_result(file_path)
        results = model.simulation.results
        self.assertEqual(len(results), 1)
        self.assertIn('Depth', results)
        self.assertEqual(results.Depth.shape, (49, 9140))


if __name__ == '__main__':
    unittest.main()

