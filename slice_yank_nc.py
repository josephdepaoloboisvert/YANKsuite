import netCDF4 as nc
import sys

def slice_yank_nc(src_nc, trg_nc, up2iter):
    src = nc.Dataset(src_nc)
    trg = nc.Dataset(trg_nc, mode='w')
    up2iter = int(up2iter)
    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        if name != 'iteration':
            trg.createDimension(name, len(dim) if not dim.isunlimited() else None)
        elif name == 'iteration':
            trg.createDimension(name, up2iter)
        else:
            raise Exception('Wut')
    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})
    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)
        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
        # Copy the variables values (as 'f4' eventually)
        trg.variables[name] = src.variables[name][:up2iter]
    # Save the file
    trg.close()

if __name__ == '__main__':
    slice_yank_nc(sys.argv[1], sys.argv[2], sys.argv[3])
    #  citation - https://coderedirect.com/questions/236272/copy-netcdf-file-using-python