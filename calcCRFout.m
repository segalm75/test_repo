% calcCRFout
% subroutine to read and calc CRF output
% files from libRadtran
%------------------------------------------
%% Details of the function:
% NAME:
%   calcCRFout
% 
% PURPOSE:
%   read relevant parameters from sfc/toa.dat files
%   calculates CRF surface/TOA and plots results
% 
%
% CALLING SEQUENCE:
%  function [d] = calcCRFout
%
% INPUT:
%  - .dat sfc/toa files
% 
% 
% OUTPUT:
%  -  d struct: CRF quantities
%  -  3D surface plots (e.g. LWP/CTK/CRF)
%
% DEPENDENCIES: 
%
% NEEDED FILES:
%  - allsky_sfc/allsky_toa, clrsky_sfc,clrsky_toa.dat
%
% EXAMPLE:
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer, NASA Ames, Nov, 3, 2014
%
% -------------------------------------------------------------------------
clear all;close all;
%% define input directory
filetmp = '*sky*.dat';
infile = getfullname__('*sky*.dat','cygwin','Select a .dat summary (sfc/toa) file');
[pname, fname, ext] = fileparts(infile);
pname = [pname, filesep];
clrfiles = dir([pname,'clr*.dat']);
cldfiles = dir([pname,'all*.dat']);
%%
%% parse file names, open files, calc CRF
for i=1:length(cldfiles)
    filestr = cldfiles(i).name;
    label   = strtok(filestr,'.');
    [dat.(label),del,headerline] = importdata([pname,filestr]);
    % match sza clear with allsky
    numsza  = unique(dat.(label).data(:,1));
    CRFtype = strsplit(label,'_');
    % calculate CRF surface
    if strcmp(CRFtype,'sfc')
            clrfile = strcat('clrsky','_',CRFtype{:},'.dat');
            dat.clr = importdata([pname,clrfile]);
            for ll=1:length(numsza)
                row = find(dat.clr.data(:,1)==numsza(ll));
                clr = dat.clr.data(row,3) - dat.clr.data(row,4);
                cld = dat.(label).data(dat.(label).data(:,1)==numsza(ll),8) - ...
                      dat.(label).data(dat.(label).data(:,1)==numsza(ll),9);
                dat.(CRFtype{:}) = clr - cld;
                % create 2D array of CRF results
                % sort lwp values
                [lwpsort, Isort] = sort(dat.(label).data(:,4),'ascend');
                ctksort = dat.(label).data(Isort,5);
                CRFsfcsort = dat.(CRFtype{:})(Isort);
                numlwp  = unique(dat.(label).data(:,4));
                numctk  = unique(dat.(label).data(:,5));
                sfc.CRF2D = reshape(CRFsfcsort,length(numlwp),length(numctk));
                figure;surf(numlwp,numctk,sfc.CRF2D' ,sfc.CRF2D') ;
                figure;scatter3(dat.(label).data(:,4),dat.(label).data(:,5),dat.(CRFtype{:}),12,dat.(CRFtype{:}));
            end
    else
    % calculate CRF TOA
            clrfile = strcat('clrsky','_',CRFtype{:},'.dat');
            dat.clr = importdata([pname,clrfile]);
            for ll=1:length(numsza)
                row = find(dat.clr.data(:,1)==numsza(ll));
                clr = dat.clr.data(row,4);
                cld = dat.(label).data(dat.(label).data(:,1)==numsza(ll),9);
                dat.(CRFtype{:}) = clr - cld;
                % create 2D array of CRF results
                numlwp  = unique(dat.(label).data(:,4));
                numctk  = unique(dat.(label).data(:,5));
                toa.CRF2D = reshape(dat.(CRFtype{:}),numlwp,numctk);
            end
    
    end
    
end

%% plot CRF 3D
% create 2D array of results (LWPxCTK)


%%












sza=[];ref=[];cbh=[];lwpind=[];lwp=[];ctk=[];lwc=[];cod=[];
sfc.dir =[];         toa.dir = [];
sfc.diffdn =[];      toa.diffdn = [];
sfc.diffup =[];      toa.diffup = [];

for i=1:length(outfiles)
    str = outfiles(i).name;
    del = '_';
    label = strsplit(str,del);
        if strcmp(filetmp,'*allsky*.out')
            sza   =[sza;   str2num(label{3}(4:end))]; 
            ref   =[ref;   str2num(label{4}(4:end))];   % [um]
            cbh   =[cbh;   str2num(label{5}(4:end))];   % [m]
            lwpind=[lwpind;str2num(label{6}(4:end))]; 
            ctk   =[ctk;   str2num(label{7}(4:end))];   % [m];
            lwp   =[lwp;   lwparray(lwpind(i))];        % [g/m2]
            lwc   =[lwc;   lwp(i)/ctk(i)];              % [g/m3]
            cod   =[cod;   (3/2)*lwp(i)/ref(i)];        % [-]
        else
            sza   =[sza;   str2num(label{3}(4:end))]; 
        end
        
        % read variables from .out file
        outtmp=load([pname,str]);
        sfc.dir    = [sfc.dir;   outtmp(1,5)];
        sfc.diffdn = [sfc.diffdn;outtmp(1,7)];
        sfc.diffup = [sfc.diffup;outtmp(1,8)];
        toa.dir    = [toa.dir;       outtmp(2,5)];
        toa.diffdn = [toa.diffdn;    outtmp(2,7)];
        toa.diffup = [toa.diffup;    outtmp(2,8)];   
end
%%
%% save to new struct and data files
%% write the surface/toa output to file/struct
if strcmp(filetmp,'*allsky*.out')
    sfcfi=[pname 'allsky_sfc']; sfcfl=[sfcfi,'.dat'];
    toafi=[pname 'allsky_toa']; toafl=[toafi,'.dat'];
    
    header=['sza' '   ' 'Reff[um]' '   ' 'CBH [m]' '   ' 'LWP [g/m2]'...
            '   ' 'CTK [m]' '   ' 'LWC [g/m3]' '   ' 'COD [-]' '   ' 'dir irrad sum [W/m2]'...
            '   ' 'diff irrad down sum [W/m2]' '   ' 'diff irrad up sum [W/m2]'];
    % surface .dat file
    disp(['writing to ascii file: ' sfcfl])
    dlmwrite(sfcfl,header,'delimiter','');
    dat=[sza,ref,cbh,lwp,ctk,lwc,sfc.dir,sfc.diffdn,sfc.diffup];
    dlmwrite(sfcfl,dat,'-append','delimiter','\t','precision',7);
    
    % surface .mat file
    disp(['saving to mat file: ' sfcfi])
    d.sza=sza;
    d.ref=ref;
    d.cbh=cbh;
    d.lwp=lwp;
    d.ctk=ctk;
    d.lwc=lwc;
    d.cod=cod;
    d.dir=sfc.dir;
    d.diffdn=sfc.diffdn;
    d.diffup=sfc.diffup;
    save(sfcfi, 'd');

    % toa file .dat
    disp(['writing to ascii file: ' toafl])
    dlmwrite(toafl,header,'delimiter','');
    dat=[sza,ref,cbh,lwp,ctk,lwc,toa.dir,toa.diffdn,toa.diffup];
    dlmwrite(toafl,dat,'-append','delimiter','\t','precision',7);
    
    % toa .mat file
    disp(['saving to mat file: ' toafi])
    d.sza=sza;
    d.ref=ref;
    d.cbh=cbh;
    d.lwp=lwp;
    d.ctk=ctk;
    d.lwc=lwc;
    d.cod=cod;
    d.dir=toa.dir;
    d.diffdn=toa.diffdn;
    d.diffup=toa.diffup;
    save(toafi, 'd');
else
    sfcfi=[pname 'clrsky_sfc']; sfcfl=[sfcfi,'.dat'];
    toafi=[pname 'clrsky_toa']; toafl=[toafi,'.dat'];
    header=['sza' '   ' 'dir irrad sum [W/m2]'...
            '   ' 'diff irrad down sum [W/m2]' '   ' 'diff irrad up sum [W/m2]'];
    % surface file
    disp(['writing to ascii file: ' sfcfl])
    dlmwrite(sfcfl,header,'delimiter','');
    dat=[sza,sfc.dir,sfc.diffdn,sfc.diffup];
    dlmwrite(sfcfl,dat,'-append','delimiter','\t','precision',7);
    
    % surface .mat file
    disp(['saving to mat file: ' sfcfi])
    d.sza=sza;
    d.dir=sfc.dir;
    d.diffdn=sfc.diffdn;
    d.diffup=sfc.diffup;
    save(sfcfi, 'd');
    
    % toa file
    disp(['writing to ascii file: ' toafl])
    dlmwrite(toafl,header,'delimiter','');
    dat=[sza,toa.dir,toa.diffdn,toa.diffup];
    dlmwrite(toafl,dat,'-append','delimiter','\t','precision',7);
    
    % toa .mat file
    disp(['saving to mat file: ' toafi])
    d.sza=sza;
    d.dir=toa.dir;
    d.diffdn=toa.diffdn;
    d.diffup=toa.diffup;
    save(toafi, 'd');
end

%%
