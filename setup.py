from setuptools import setup, find_packages

setup(
    name='AdhModel',
    version='0.4.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'netcdf4',
        'jupyter',
        'param>=1.9.0',
        'holoviews>=1.12.3',
        'datashader>=0.7.0',
        'geoviews>=1.6.2',
        'panel>=0.6.0',
        'bokeh>=1.1.0',
        'cartopy>=0.17.0',
        'xarray>=0.11.0',
        'colorcet>=1.0.0',
        'notebook>=5.5.0',
        'fiona',
        'gdal>=2.4.1',
        'rasterio>=1.0.24',
        'xmscore',
        'xmsinterp',
        'xmsgrid',
        'xmsmesh',
        'opencv',
        'quest>=3.1.1',
        'ulmo>=0.8.5',
        'lancet>=0.9.0',
        'pyct >=0.4.4',
        'setuptools >=30.3.0',
        'pyviz_comms >=0.6.0',
        'nodejs >=9.11.1',
        'earthsim>=1.1.1a2',
        'jupyterlab',
        'pytest',
        'pyflakes',
        'nbsmoke',
        'genesis>=0.0.5'
    ],
    description="Parameterized class for reading/write/specifying all aspects of an Adaptive Hydraulics Model",
    url="https://github.com/erdc/AdhModel",
    )
