% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023

function [C, varargout] = bmCoilSense_cartesian_secondary(y, C, y_ref, C_ref, N_u, dK_u, convCond)

% initial -----------------------------------------------------------------

nIterSmooth = 2; % ----------------------------------------------------------- magic number

if not(strcmp(class(y), 'single'))
    error('y must be of class ''single'' .');
    return; 
end
if not(strcmp(class(y_ref), 'single'))
    error('y must be of class ''single'' .');
    return; 
end

N_u     = double(N_u(:)');
dK_u    = double(dK_u(:)');
nPt     = double(prod(N_u(:))); 
imDim   = size(N_u(:), 1);
ve      = prod(dK_u(:))*ones(1, nPt);  

nCh_array   = size(y, 2); 
C           = C/nCh_array; 
y_ref       = nCh_array*y_ref/bmY_norm(y_ref, ve)*mean(bmCol(bmY_norm(y, ve, false))); 

y           = cat(2, y_ref(:), y);   
C           = cat(2, C_ref(:), bmColReshape(C, N_u)); 
x           = bmColReshape(bmFourier_inv(y, N_u, dK_u, C),  N_u); 
F_inv_y     = bmColReshape(bmFourier_inv(y, N_u, dK_u),     N_u); 
nCh         = size(y, 2); 
ve          = bmY_ve_reshape(ve, size(y)); 
myZero      = zeros(prod(N_u(:)), 1); 

% END_initial -------------------------------------------------------------


while convCond.check()
    
    v = bmFourier(x, N_u, dK_u, C) - y;
    w = bmFourier_conjTrans(ve.*v, N_u, dK_u);
    d_C = 2*repmat(conj(x), [1, nCh]).*w;
    d_C(:, 1) = myZero; 
    Ad_C = bmFourier(x, N_u, dK_u, d_C);
    lambda_C = real(  Ad_C(:)' * (ve(:).*v(:))  )/real(  Ad_C(:)' * (ve(:).*Ad_C(:))  );
    
    x = bmColReshape(bmCoilSense_pinv(C, F_inv_y, N_u), N_u); 
    C = C - lambda_C*d_C;
      
    convCond.disp_info('bmCoilSense_cartesian_secondary'); 
end

% final -------------------------------------------------------------------

C = C(:, 2:end); 
nCh = size(C, 2); 
C = bmBlockReshape(C, N_u);
x = bmBlockReshape(x, N_u); 

for i = 1:nCh
    if imDim == 1
        C(:, i) = bmImPseudoDiffusion(C(:, i), nIterSmooth);
    elseif imDim == 2
        C(:, :, i) = bmImPseudoDiffusion(C(:, :, i), nIterSmooth);
    elseif imDim == 3
        C(:, :, :, i) = bmImPseudoDiffusion(C(:, :, :, i), nIterSmooth);
    end
end

C = C*nCh_array; 

if nargout > 1
    varargout{1} = convCond; 
end
if nargout > 2
    varargout{2} = x; 
end
% END_final ---------------------------------------------------------------

end