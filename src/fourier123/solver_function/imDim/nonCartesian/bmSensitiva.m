% Bastien Milani
% December 2025
% Martigny - Switzerland

function x = bmSensitiva(   x, y, ve, C, Gu, Gut, frSize, ...
                            delta, nCGD, ve_max, ...
                            regularizer, ...
                            nIter, witnessInfo)

% initial -----------------------------------------------------------------
myEps       = 10*eps('single'); % --------------------------------------------- magic_number

y           = single(y);  
nCh         = size(y, 2); 
N_u         = double(single(Gu.N_u(:)')); 
frSize      = double(single(frSize(:)'));

dK_u        = double(single(Gu.d_u(:)'));
C           = single(bmColReshape(C, N_u));  
ve          = single(bmY_ve_reshape(ve, size(y)));  

KFC         = single(bmColReshape(bmKF(          C,  N_u, frSize, dK_u, nCh, Gu.kernel_type, Gu.nWin, Gu.kernelParam), frSize)); 
KFC_conj    = single(bmColReshape(bmKF_conj(conj(C), N_u, frSize, dK_u, nCh, Gu.kernel_type, Gu.nWin, Gu.kernelParam), frSize)); 

x           = single(bmColReshape(x, frSize)); 

dX_u        = single(  (1./single(dK_u))./single(N_u)  );
deltaR      = single(  prod(dX_u(:))); 
HX          = single(  prod(dX_u(:))); 
HZ          = single(  prod(dX_u(:))); 


if isempty(ve_max)
   ve_max = max(ve(:));  
end
HY = min(ve, single(ve_max)); 

private_init_witnessInfo(witnessInfo, nIter, 'sensa', frSize, N_u, dK_u, nCGD, ve_max); 


% monitoring ----------------------------------------------------------
c = 1; 

temp_r                      = y - private_M(x, Gu, frSize, KFC); 
dafi                        = real(  temp_r(:)'*(HY(:).*temp_r(:))  );

temp_r                      = private_F(x, regularizer, frSize, dX_u);
regul                       = deltaR*real(  temp_r(:)'*temp_r(:)  );

witnessInfo.param{7}(1, 1)  = dafi;
witnessInfo.param{8}(1, 1)  = regul;
witnessInfo.param{9}(1, 1)  = 0.5*dafi + 0.5*regul;

witnessInfo.watch(c, x, frSize, 'initial');
% END_monitoring ------------------------------------------------------


% END_initial -------------------------------------------------------------

% main_loop ---------------------------------------------------------------
for c = 1:nIter
    
    % L'Aube
    res_y_next          = y - private_M(x, Gu, frSize, KFC);
    res_z_next          = 0 - private_F(x, regularizer, frSize, dX_u); 
    dagM_res_y_next     = private_dagM(res_y_next, Gut, HX, HY, frSize, KFC_conj);
    dagF_res_z_next     = private_dagF(res_z_next, regularizer, HX, HZ, frSize, dX_u); 
    dagA_res_next       = dagM_res_y_next + delta*dagF_res_z_next;  
    p_next              = dagA_res_next;
    sqn_dagA_res_next   = real(   dagA_res_next(:)'*(HX.*dagA_res_next(:))   );
    
    
    for i = 1:nCGD
        
        % Le Matin
        res_y_curr          = res_y_next;
        res_z_curr          = res_z_next;
        sqn_dagA_res_curr   = sqn_dagA_res_next; 
        p_curr              = p_next;
        
        if (sqn_dagA_res_curr < myEps) 
            break;
        end
        
        % Le Midi
        Mp_curr             = private_M(p_curr, Gu, frSize, KFC);
        Fp_curr             = private_F(p_curr, regularizer, frSize, dX_u);
        sqn_Mp_curr         = real(   Mp_curr(:)'*(HY(:).*Mp_curr(:))   );
        sqn_Fp_curr         = real(   Fp_curr(:)'*(delta*HZ(:).*Fp_curr(:))   );
        sqn_Ap_curr         = sqn_Mp_curr + sqn_Fp_curr; 

        % Le Soir
        a   = sqn_dagA_res_curr/sqn_Ap_curr;
        x = x + a*p_curr;
                
        if (i == nCGD)
           break;  
        end
        
        % La Nouvelle Aube
        res_y_next          = res_y_curr - a*Mp_curr;
        res_z_next          = res_z_curr - a*Fp_curr;

        dagM_res_y_next       = private_dagM(res_y_next, Gut, HX, HY, frSize, KFC_conj);
        dagF_res_z_next       = private_dagF(res_z_next, regularizer, HX, HZ, frSize, dX_u);
        
        dagA_res_next       = dagM_res_y_next + delta*dagF_res_z_next; 
        sqn_dagA_res_next   = real(   dagA_res_next(:)'*(HX.*dagA_res_next(:))   );
        
        b                   = sqn_dagA_res_next/sqn_dagA_res_curr; 
        
        p_next              = dagA_res_next + b*p_curr;
        
    end % end CGD
    
    % monitoring ----------------------------------------------------------
    temp_r                          = y - private_M(x, Gu, frSize, KFC);
    dafi                            = real(  temp_r(:)'*(HY(:).*temp_r(:))  );

    temp_r                          = private_F(x, regularizer, frSize, dX_u);
    regul                           = deltaR*real(  temp_r(:)'*temp_r(:)  );

    witnessInfo.param{7}(1, c+1)    = dafi; 
    witnessInfo.param{8}(1, c+1)    = regul; 
    witnessInfo.param{9}(1, c+1)    = 0.5*dafi + 0.5*regul;
    witnessInfo.watch(c, x, frSize, 'loop'); 
    % END_monitoring ------------------------------------------------------

end
% END_main_loop -----------------------------------------------------------

% final -------------------------------------------------------------------
witnessInfo.watch(c, x, frSize, 'final'); 
x = bmBlockReshape(x, frSize); 
% END_final ---------------------------------------------------------------

end



function private_init_witnessInfo(witnessInfo, nIter, argName, frSize, N_u, dK_u, nCGD, ve_max)

witnessInfo.param_name{1}    = 'recon_name'; 
witnessInfo.param{1}         = argName; 

witnessInfo.param_name{2}    = 'dK_u'; 
witnessInfo.param{2}         = dK_u; 

witnessInfo.param_name{3}    = 'N_u'; 
witnessInfo.param{3}         = N_u; 

witnessInfo.param_name{4}    = 'frSize'; 
witnessInfo.param{4}         = frSize; 

witnessInfo.param_name{5}    = 'nCGD'; 
witnessInfo.param{5}         = nCGD; 

witnessInfo.param_name{6}    = 've_max'; 
witnessInfo.param{6}         = ve_max;

witnessInfo.param_name{7}    = 'residu'; 
witnessInfo.param{7}         = zeros(1, nIter);

witnessInfo.param_name{8}    = 'regul'; 
witnessInfo.param{8}         = zeros(1, nIter);

witnessInfo.param_name{9}    = 'objective'; 
witnessInfo.param{9}         = zeros(1, nIter);

end


% model_and_regularizer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% forward_model
function M_x = private_M(x, Gu, frSize, KFC)
    M_x     = bmShanna(x, Gu, KFC, frSize, 'MATLAB');
end


% adjoint_model
function dagM_y = private_dagM(y, Gut, HX, HY, frSize, KFC_conj)
    dagM_y = (1./HX).*bmNakatsha(HY.*y, Gut, KFC_conj, true, frSize, 'MATLAB'); % negative_gradient
end


% forward_regularizer
function F_x = private_F(x, regularizer, frSize, dX_u)
    if strcmp(regularizer, 'spatial_derivative')
        F_x = bmBackGradient(x, frSize, dX_u); 
    elseif strcmp(regularizer, 'identity')
        F_x = x;
    end
end


% adjoint_regularizer
function dagF_z = private_dagF(z, regularizer, HX, HZ, frSize, dX_u)

    if strcmp(regularizer, 'spatial_derivative')
        dagF_z = (1./HX).*bmBackGradientT(HZ.*z, frSize, dX_u); 
    elseif strcmp(regularizer, 'identity')
        dagF_z = z; 
    end

end
% END_model_and_regularizer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




