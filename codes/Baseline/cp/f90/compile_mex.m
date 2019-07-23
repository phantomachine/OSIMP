ArrayFlag = '-compatibleArrayDims';

eval(['mex -c mod_parameters.f90 ',ArrayFlag]);
eval(['mex -c mod_spline.f90 ',ArrayFlag]);
eval(['mex -c mod_foccp.f90 ',ArrayFlag]);
eval(['mex -c mod_calccp.f90 ',ArrayFlag]);

if (ispc)
    eval(['mex calccp.f90 mod_calccp.obj mod_foccp.obj mod_spline.obj mod_parameters.obj ',ArrayFlag]);
else
    eval(['mex calccp.f90 mod_calccp.o mod_foccp.o mod_spline.o mod_parameters.o ',ArrayFlag]);
end
