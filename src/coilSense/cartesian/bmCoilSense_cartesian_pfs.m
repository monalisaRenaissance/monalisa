% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023

function y_pfs = bmCoilSense_cartesian_pfs(y, im_ref, N_u, dK_u, myPerm, myFlip, myShift)

y               = bmColReshape(y, N_u); 
N_u             = N_u(:)'; 
dK_u            = dK_u(:)';
imDim           = size(N_u(:), 1); 
myPerm          = myPerm(:)';
myFlip          = myFlip(:)';
myShift         = myShift(:)'; 

x = bmBlockReshape(bmFourier_inv(y, N_u, dK_u), N_u); 
x = permute(x, [myPerm, imDim+1]); 
for i = 1:imDim
   if myFlip(1, i) > 0
        x = flip(x, i); 
   end
end
x = circshift(x, myShift); 

y_pfs = bmColReshape(bmFourier(x, N_u, dK_u), N_u); 


x_rms   = bmRMS(x, N_u); 
x_rms   = bmBlockReshape(x_rms, N_u); 
myImage = cat(imDim + 1, x_rms, im_ref);
bmImage(myImage/mean(abs(myImage(:)))); 

end