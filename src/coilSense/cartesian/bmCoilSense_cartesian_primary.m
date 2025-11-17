% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023

function C = bmCoilSense_cartesian_primary(y, y_ref, C_ref, N_u, dK_u, m)

nIter_smooth = 2; % ----------------------------------------------------------- magic number
L_nIter = 1000; % ------------------------------------------------------------- magic number
L_th = 1e-4; % ---------------------------------------------------------------- magic number

N_u         = double(N_u(:)');
dK_u        = double(dK_u(:)'); 
nPt         = double(prod(N_u(:)));
imDim       = size(N_u(:), 1);  
nCh         = size(y(:), 1)/nPt;
ve          = prod(dK_u(:))*ones(1, nPt); 

nCh_array   = size(y, 2);  
y_ref       = nCh_array*y_ref/bmY_norm(y_ref, ve)*mean(bmCol(bmY_norm(y, ve, false))); 

x_ref       = bmBlockReshape(bmFourier_inv(y_ref, N_u, dK_u), N_u);
x           = bmBlockReshape(bmFourier_inv(y, N_u, dK_u),     N_u); 
 
m           = logical(bmBlockReshape(m, N_u)); 
m_neg       = not(m); 

anat_ref    = x_ref./C_ref; 
anat_ref(m_neg) = 1; 

z = zeros([N_u, nCh], 'single'); 
C = complex(z, z); 

for i = 1:nCh
    if imDim == 1
        temp_im = x(:, i); 
        temp_im = bmImPseudoDiffusion_inMask(temp_im./anat_ref, m, nIter_smooth); 
        temp_im(m_neg) = 0; 
        C(:, i) = bmImLaplaceEquationSolver(temp_im, m, L_nIter, L_th, 'omp');
        
    elseif imDim == 2
        temp_im = x(:, :, i);  
        temp_im = bmImPseudoDiffusion_inMask(temp_im./anat_ref, m, nIter_smooth);
        temp_im(m_neg) = 0; 
        C(:, :, i) = bmImLaplaceEquationSolver(temp_im, m, L_nIter, L_th, 'omp');

    elseif imDim == 3
        temp_im = x(:, :, :, i);  
        temp_im = bmImPseudoDiffusion_inMask(temp_im./anat_ref, m, nIter_smooth);
        temp_im(m_neg) = 0; 
        C(:, :, :, i) = bmImLaplaceEquationSolver(temp_im, m, L_nIter, L_th, 'omp');
        
    end
end

C = C*nCh_array; 
C = bmBlockReshape(C, N_u); 


end