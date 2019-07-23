ArrayFlag = '-compatibleArrayDims';

eval(['mex -c mod_parameters.f90 ',ArrayFlag]);
eval(['mex -c mod_spline.f90 ',ArrayFlag]);
eval(['mex -c mod_focpar.f90 ',ArrayFlag]);
eval(['mex -c mod_calcpar.f90 ',ArrayFlag]);

if (ispc)
    eval(['mex calcpar.f90 mod_calcpar.obj mod_focpar.obj mod_spline.obj mod_parameters.obj ',ArrayFlag]);
else
    eval(['mex calcpar.f90 mod_calcpar.o mod_focpar.o mod_spline.o mod_parameters.o ',ArrayFlag]);
end
