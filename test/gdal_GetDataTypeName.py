from osgeo import gdal



for i in range(15):
    dtype = gdal.GetDataTypeName(i)
    print(f'{i}: {dtype}')