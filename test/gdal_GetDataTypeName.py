from osgeo import gdal


def get_gdal_data_type(label):
    data_type_mapping = {
        'Unknown': gdal.GDT_Unknown,
        'Byte': gdal.GDT_Byte,
        'UInt16': gdal.GDT_UInt16,
        'Int16': gdal.GDT_Int16,
        'UInt32': gdal.GDT_UInt32,
        'Int32': gdal.GDT_Int32,
        'Float32': gdal.GDT_Float32,
        'Float64': gdal.GDT_Float64,
        'CInt16': gdal.GDT_CInt16,
        'CInt32': gdal.GDT_CInt32,
        'CFloat32': gdal.GDT_CFloat32,
        'CFloat64': gdal.GDT_CFloat64,
        'UInt64': gdal.GDT_UInt64,
        'Int64': gdal.GDT_Int64,
        'Int8': gdal.GDT_Int8,
    }

    if label in data_type_mapping:
        return data_type_mapping[label]
    else:
        raise ValueError("Invalid GDAL data type label")


for i in range(15):
    dtype = gdal.GetDataTypeName(i)
    print(get_gdal_data_type(dtype))
    # print(f'{i}: {dtype}')