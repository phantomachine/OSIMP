ArrayFlag = '-compatibleArrayDims';

eval(['mex -c mod_parameters.f90 ',ArrayFlag]);
eval(['mex -c mod_spline.f90 ',ArrayFlag]);
eval(['mex -c mod_foc.f90 ',ArrayFlag]);
eval(['mex -c mod_main.f90 ',ArrayFlag]);

eval(['mex -c mod_focnc.f90 ',ArrayFlag]);
eval(['mex -c mod_calcnc.f90 ',ArrayFlag]);

if (ispc)
    eval(['mex sus_iter.f90 mod_main.obj mod_foc.obj mod_spline.obj mod_parameters.obj ',ArrayFlag]);
    eval(['mex calcnc.f90 mod_calcnc.obj mod_focnc.obj mod_spline.obj mod_parameters.obj ',ArrayFlag]);
else
    eval(['mex calcnc.f90 mod_calcnc.o mod_focnc.o mod_spline.o mod_parameters.o ',ArrayFlag]);
    eval(['mex sus_iter.f90 mod_main.o mod_foc.o mod_spline.o mod_parameters.o ',ArrayFlag]);
end
