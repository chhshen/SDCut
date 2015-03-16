function I = sp_load_image(image_fname)

I = imread(image_fname);
if ndims(I) == 3
    I = double(ow_rgb2gray(I))./255;
else
    I = double(I)./255;
end
