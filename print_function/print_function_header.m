function fid=print_function_header(fcn_name,out,in,comment)
%PRINT_FUNCTION_HEADER

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%11/30/03 - updated for the CDC
%Peter Minh 2016 - update for optimization
comment=['%   [',upper(out),'] = ',...
	 upper(fcn_name),'(',upper(in),')  ',comment];
    
comment = fixlength(comment,' ',150,'   ');

%disp(['[creating ',upper(fcn_name),'.m]']);
fid = fopen(['m_files/',fcn_name,'.m'],'w');
fprintf(fid,['function [',out,'] = %s(',in,')\n'],fcn_name);
fprintf(fid,['%% %s\n'],upper(fcn_name));
fprintf(fid,['%s\n\n'],comment);
fprintf(fid,'%%Eric Westervelt\n');
fprintf(fid,'%%2016 Version: Peter Minh\n');
fprintf(fid,'%%%s\n\n',datestr(now));
