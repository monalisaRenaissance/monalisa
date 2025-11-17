% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023

function x = bmSteva(  x, z, u, y, ve, C, Gu, Gut, frSize, ...
                               delta, rho, nCGD, ve_max, ...
                               nIter, witnessInfo)

% initial -----------------------------------------------------------------
myEps   = 10*eps('single'); % -------------------------------------------------- magic number
    
y               = single(y);   
nCh             = size(y, 2);
N_u             = double(int32(Gu.N_u(:)'));
frSize          = double(frSize(:)');
if isempty(frSize)
   frSize = N_u;  
end
nPt_u           = prod(frSize(:)); 
imDim           = size(N_u(:), 1);  
dK_u            = double(single(Gu.d_u(:)'));
dX_u            = single(  (1./single(dK_u))./single(N_u)  );

if isempty(ve_max)
   ve_max = max(ve(:));  
end


C               = single(bmColReshape(C, N_u));
KFC             = single(bmKF(          C,  N_u, frSize, dK_u, nCh, Gu.kernel_type, Gu.nWin, Gu.kernelParam));
KFC_conj        = single(bmKF_conj(conj(C), N_u, frSize, dK_u, nCh, Gu.kernel_type, Gu.nWin, Gu.kernelParam));

[delta, rho]    = private_init_delta_rho(delta, rho, nIter); 

% HX            = single(  prod(dX_u(:)));
% HZ            = HX; 
HX              = single(  prod(dX_u(:))*mean(abs(bmColReshape(C, frSize)).^2, 2)  ); % Image-space preconditioning of K. Pruessman. 
HZ              = single(  repmat(HX(:), [1, imDim])  );
HY              = min(single(  bmY_ve_reshape(ve, size(y))  ), single(ve_max)); 
deltaR          = single(  prod(dX_u(:)));

x = single(bmColReshape(x, frSize));

if isempty(z)
    z = private_F(x, frSize, dX_u); 
end
if isempty(u)
    u = bmZero([nPt_u, imDim], 'complex_single');
end

private_init_witnessInfo(witnessInfo, nIter, 'steva', frSize, N_u, dK_u, delta, rho, nCGD, ve_max); 

% monitoring --------------------------------------------------------------
    temp_r          = y - private_M(x, Gu, frSize, KFC);
    dafi            = real(  temp_r(:)'*(HY(:).*temp_r(:))  ); 
    
    temp_r          = private_F(x, frSize, dX_u);
    regul           = deltaR*sum(abs(  real(temp_r(:))  )) + deltaR*sum(abs(  imag(temp_r(:))  ));
    
    c               = 1;
    
    objective = 0.5*dafi + 0.5*delta(1, c)*regul;
    
    witnessInfo.param{9}(1, c)      = dafi; 
    witnessInfo.param{10}(1, c)     = regul; 
    witnessInfo.param{11}(1, c)     = objective;
    witnessInfo.watch(c, x, frSize, 'initial');
% END_monitoring ----------------------------------------------------------


% END_initial -------------------------------------------------------------



% ADMM loop ---------------------------------------------------------------
for c = 1:nIter

    res_y_next   = y - private_M(x, Gu, frSize, KFC);
    res_z_next   = (z - u) - private_F(x, frSize, dX_u);
    
    dagM_res_y_next  =  private_dagM(res_y_next, Gut, HX, HY, frSize, KFC_conj); 
    dagF_res_z_next  = rho(1, c)*private_dagF(res_z_next, HX, HZ, frSize, dX_u);
                                
    
    dagA_res_next   = dagM_res_y_next + dagF_res_z_next; 
    p_next          = dagA_res_next; 
    
    sqn_dagA_res_next = real(   dagA_res_next(:)'*(HX.*dagA_res_next(:))   );
    

    for i = 1:nCGD
        
        res_y_curr   = res_y_next;
        res_z_curr   = res_z_next;
        sqn_dagA_res_curr = sqn_dagA_res_next; 
        p_curr = p_next; 
        
        
        if (sqn_dagA_res_curr < myEps)
            break;
        end
        
        
        Mp_curr     = private_M(p_curr, Gu, frSize, KFC);
        Fp_curr     = private_F(p_curr, frSize, dX_u);
        
        sqn_Mp_curr      = real(   Mp_curr(:)'*(HY(:).*Mp_curr(:))   );
        sqn_Fp_curr      = real(   Fp_curr(:)'*(rho(1, c)*HZ(:).*Fp_curr(:))   );
        sqn_Ap_curr       = sqn_Mp_curr + sqn_Fp_curr;
        
        a   = sqn_dagA_res_curr/sqn_Ap_curr;
        
        x = x + a*p_curr;
        
        if i == nCGD
            break;
        end
        
        res_y_next          = res_y_curr - a*Mp_curr;
        res_z_next          = res_z_curr - a*Fp_curr;
        
        dagM_res_y_next   = private_dagM(res_y_next, Gut, HX, HY, frSize, KFC_conj);
        dagF_res_z_next   = rho(1, c)*private_dagF(res_z_next, HX, HZ, frSize, dX_u); 
        
        dagA_res_next     = dagM_res_y_next + dagF_res_z_next;
        sqn_dagA_res_next = real(   dagA_res_next(:)'*(HX.*dagA_res_next(:))   );
        
        b = sqn_dagA_res_next/sqn_dagA_res_curr; 
        
        p_next           = dagA_res_next + b*p_curr;

    end

    bGx_plus_u      = private_F(x, frSize, dX_u) + u; 
    z               = bmProx_oneNorm(bGx_plus_u, delta(1, c)*deltaR/rho(1, c)./HZ(:)); 
    u               = bGx_plus_u - z; 
    
    
    % monitoring ----------------------------------------------------------
    temp_r          = y - private_M(x, Gu, frSize, KFC);
    dafi            = real(  temp_r(:)'*(HY(:).*temp_r(:))  ); 
    
    temp_r          = private_F(x, frSize, dX_u);
    regul           = deltaR*sum(abs(  real(temp_r(:))  )) + deltaR*sum(abs(  imag(temp_r(:))  ));
    
    objective       = 0.5*dafi + 0.5*delta(1, c)*regul;
    
    witnessInfo.param{9}(1, c+1)  = dafi; 
    witnessInfo.param{10}(1, c+1) = regul; 
    witnessInfo.param{11}(1, c+1) = objective; 

    witnessInfo.watch(c, x, frSize, 'loop'); 
    % END_monitoring ------------------------------------------------------
    
end
% END_ADMM loop -----------------------------------------------------------

% final -------------------------------------------------------------------
witnessInfo.watch(c, x, frSize, 'final'); 
x = bmBlockReshape(x, frSize);
% END_final ---------------------------------------------------------------

end

function [delta, rho] = private_init_delta_rho(delta, rho, nIter)

rho             = single(  abs(rho(:))  );
delta           = single(  abs(delta(:))  );
if size(delta, 1) == 1
    delta   = linspace(delta, delta, nIter);
elseif size(delta, 1) == 2
    delta   = linspace(delta(1, 1), delta(2, 1), nIter);
end
delta = delta(:)';
if size(rho, 1) == 1
    rho     = linspace(rho,     rho, nIter);
elseif size(rho, 1) == 2
    rho     = linspace(rho(1, 1),     rho(2, 1), nIter);
end
rho = rho(:)';

end

function private_init_witnessInfo(witnessInfo, nIter, argName, frSize, N_u, dK_u, delta, rho, nCGD, ve_max)

witnessInfo.param_name{1}    = 'recon_name'; 
witnessInfo.param{1}         = argName; 

witnessInfo.param_name{2}    = 'dK_u'; 
witnessInfo.param{2}         = dK_u; 

witnessInfo.param_name{3}    = 'N_u'; 
witnessInfo.param{3}         = N_u; 

witnessInfo.param_name{4}    = 'frSize'; 
witnessInfo.param{4}         = frSize; 

witnessInfo.param_name{5}    = 'delta'; 
witnessInfo.param{5}         = delta; 

witnessInfo.param_name{6}    = 'rho'; 
witnessInfo.param{6}         = rho; 

witnessInfo.param_name{7}    = 'nCGD'; 
witnessInfo.param{7}         = nCGD; 

witnessInfo.param_name{8}    = 've_max'; 
witnessInfo.param{8}         = ve_max;

witnessInfo.param_name{9}    = 'residu'; 
witnessInfo.param{9}         = zeros(1, nIter); 

witnessInfo.param_name{10}   = 'total_variation'; 
witnessInfo.param{10}        = zeros(1, nIter); 

witnessInfo.param_name{11}   = 'objective'; 
witnessInfo.param{11}        = zeros(1, nIter); 

end




% model_and_sparsifier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% forward_model
function M_x = private_M(x, Gu, frSize, KFC)
    M_x     = bmShanna(x, Gu, KFC, frSize, 'MATLAB');
end

% adjoint_model
function dagM_y = private_dagM(y, Gut, HX, HY, frSize, KFC_conj)
    dagM_y = (1./HX).*bmNakatsha(HY.*y, Gut, KFC_conj, true, frSize, 'MATLAB'); % negative_gradient
end



% forward_sparsifier
function F_x = private_F(x, frSize, dX_u)
    F_x     = bmBackGradient(x, frSize, dX_u); 
end


% adjoint_sparsifier
function dagF_z = private_dagF(z, HX, HZ, frSize, dX_u)
    dagF_z = (1./HX).*bmBackGradientT(HZ.*z, frSize, dX_u); 
end


% END_model_and_sparsifier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




