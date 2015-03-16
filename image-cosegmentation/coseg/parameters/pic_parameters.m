
global npics

% NUMBER OF PICS
if ~exist('n_pics')
    fprintf('ERROR unspecified number of images --- set nb of im = 2\n')
    n_pics=2;
elseif ischar(n_pics)==1
        n_pics=str2num(n_pics);
end
npics=n_pics;
param.pic.pics=n_pics;


global ndecors
if ~exist('n_decor')
    n_decor=0;
elseif ischar(n_decor)==1
    n_decor=str2num(n_decor);
end
ndecors=n_decor;
param.pic.decor=n_decor;
