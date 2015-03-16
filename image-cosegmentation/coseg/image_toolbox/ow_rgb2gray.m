function res = ow_rgb2gray(im)
imd = double(im);
if size(im,3) > 1
%    0.298936021293776 % values from rgb2gray
%    0.587043074451121
%    0.114020904255103
   res = cast(0.298936021293776*imd(:,:,1) +...
         0.587043074451121*imd(:,:,2) + 0.114020904255103*imd(:,:,3),...
        class(im));
else
   res = im;
end
if strcmp(class(im),'double')
   res = min(max(res,0),1);
end