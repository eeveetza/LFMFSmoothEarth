% This script is used to validate the MATLAB/Octave implementation of
% LFMF-SmoothEarth against the values obtained using the reference C++
% implementation from NTIA on a set of different tests (available in the
% file ValidationExampleLFMFSmoothEarth.csv. As the data in this file is
% written with precision of %.6f, the test passes if the tolerance is not
% greater than 1e-6.

% Rev   Date        Author                          Description
%-------------------------------------------------------------------------------
% v1   4JUN21     Ivica Stevanovic, OFCOM         Initial implementation


clear all
close all
fclose all;
clc


% immediate printing to command window in octave
% if (isOctave)
%     page_screen_output(0);
%     page_output_immediately(1);
% end


%% begin code

% tolerance 
tol = 1e-6;

filename= 'ValidationExampleLFMFSmoothEarth.csv';


fid=fopen(filename,'r');
if (fid==-1)
    
    error('File cannot be found.');
    
end

% First line is header
readLine = fgetl(fid);
cnt_pass = 0;
cnt_fail = 0;

    
fprintf(1,'%20s %20s %20s \n', 'MATLAB', 'REF', 'DELTA' );

while(1)
    readLine = fgetl(fid);
    if (~ischar(readLine))
        break
    end
    
    
    dummy=regexp(readLine,',','split');
    htx = str2double(dummy(1));
    hrx = str2double(dummy(2));
    fmhz = str2double(dummy(3));
    ptx = str2double(dummy(4));
    ns = str2double(dummy(5));
    dkm = str2double(dummy(6));
    eps = str2double(dummy(7));
    sigma = str2double(dummy(8));
    pol = str2double(dummy(9));
    
    tl_ref = str2double(dummy(10));
    result = tl_LFMFSmoothEarth(htx, hrx, fmhz, ptx, ns, dkm, eps, sigma, pol);
    tl = result.A_btl__db;
    
    delta = abs(tl-tl_ref);
    
    
    if (abs(delta)>tol)
        cnt_fail = cnt_fail + 1;
    else
        cnt_pass = cnt_pass + 1;
    end
    
    fprintf(1,'%20.6f %20.6f %20.6f \n', tl, tl_ref, delta );
    
    
end
fclose(fid);

fprintf(1,'Successfully passed %d out of %d tests\n', cnt_pass, cnt_pass+cnt_fail);
