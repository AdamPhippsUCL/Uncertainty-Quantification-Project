% Script to define and make dMRI scheme structure

clear;
projectfolder = pwd;
schemesfolder = fullfile(projectfolder, 'Simulation Feasibility Study', 'Schemes');

%% Define sequences

% V = [delta, DELTA, bval, TE, N0, Nb]
% 
% Where: N0 - number of b=0 averages, Nb - number of b>0 averages
% (additional factor)

% Classic
b3000 = [20, 43.4, 3000, 87, 6, 3*3];
b2000 = [14, 37.4, 2000, 75, 6, 3*3];
b1500 = [23.5, 46.9, 1500, 94, 6, 3*3 ];
b500 = [10.5, 33.9, 500, 68, 6, 3*2];
b90 = [3.5, 26.9, 90, 54, 4, 3*1];

% % Fast at Classic Geom
% b50_ClassicGeom = [11.0, 34.4, 50];
% b1000_ClassicGeom = [11.0, 34.4, 1000];
% b1800_ClassicGeom = [18.0, 41.4, 1800];
% 
% Fast
b1000_Fast = [21.0, 34.9, 1000, 69, 6, 3*3 ];
b1800_Fast = [24.0, 46.3, 1800, 83, 6, 3*3 ];



%% Build schemes

schemenames = {...
    'Fast',...
};


for sindx = 1:length(schemenames)

    schemename = schemenames{sindx};
    
    switch schemename
    
        case 'Classic'
            
            Vs = [b3000; b2000; b1500; b500; b90];

        case 'Fast'

            Vs = [b1800_Fast; b1000_Fast];
    
    end
    
    scheme = BuildScheme(Vs, schemename);
    
    save(fullfile(schemesfolder, [schemename '.mat']), "scheme");

end


%% Function

function scheme = BuildScheme(scans, schemename)

for scanIndx = 1:size(scans,1)

    delta = scans(scanIndx, 1);
    Delta = scans(scanIndx,2);
    bval = scans(scanIndx,3);
    TE = scans(scanIndx,4);
    N0 = scans(scanIndx,5);
    Nb = scans(scanIndx,6);

    scheme(scanIndx).delta = delta;
    scheme(scanIndx).DELTA = Delta ;
    scheme(scanIndx).bval = bval;
    scheme(scanIndx).G = stejskal(delta,Delta,bval=bval);
    scheme(scanIndx).TE = TE;
    scheme(scanIndx).N0 = N0;
    scheme(scanIndx).Nb = Nb;
    
    scheme(scanIndx).schemename = schemename;

end
    
end


