function [one_pic_aj]=find_indices_image(lWaj,D)
% function [one_pic_aj]=find_indices_image(lWaj,D)
% create a NxI matrix M where N is the number of pixels and I the number of
% images, with M_{ij}=1 if a pixel i is in iamge j and other otherwise

I=(1:sum(lWaj))';
J=0.*I;
V=J;

lWtot                                       =1;

for iaj=1:length(lWaj)  
   J(lWtot:lWtot+lWaj(iaj)-1) = iaj;
   V(lWtot:lWtot+lWaj(iaj)-1) = D(lWtot:lWtot-1+lWaj(iaj),1);
   lWtot                                    =lWtot+lWaj(iaj);
end

one_pic_aj =sparse(I,J,V);
clear I J V
