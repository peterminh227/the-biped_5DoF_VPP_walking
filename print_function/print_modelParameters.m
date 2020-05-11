function print_modelParameters(fid)
%PRINT_MODELPARAMETERS

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%11/30/03 - updated for the CDC

fprintf(fid,'modparams=modelParameters;\n');
fprintf(fid,'g=modparams.g;\n');
fprintf(fid,'L1=modparams.L1;   L3=modparams.L3;   L4=modparams.L4;\n');
fprintf(fid,'M1=modparams.M1;   M3=modparams.M3;   M4=modparams.M4;\n');
fprintf(fid,'MY1=modparams.MY1; MY3=modparams.MY3; MY4=modparams.MY4;\n');
fprintf(fid,'MZ1=modparams.MZ1; MZ3=modparams.MZ3; MZ4=modparams.MZ4;\n');
fprintf(fid,'XX1=modparams.XX1; XX3=modparams.XX3; XX4=modparams.XX4;\n');

