% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023



function relax_map = bmIR_fit(argImagesTable, TIR, TR, varargin)


% initial -----------------------------------------------------------------
argSize         = size(argImagesTable);
s1              = prod(argSize(1:end-1));
s2              = argSize(1, end);

imagesTable     = reshape(argImagesTable, [s1, s2]);
TIR             = TIR(:)';

amp_map         = zeros(s1, 1);
relax_map       = zeros(s1, 1);


if not(size(TIR(:), 1) == s2)
    T1 = [];
    errordlg('Wrong list of arguments');
    return;
end

errorTh             = [];
lowerBound          = [];
upperBound          = [];
lsqLowerBound       = [];
lsqUpperBound       = [];
lsqcurvefitFlag     = 0;
display_flag        = false; 

if length(varargin) == 0
    1+1;
elseif length(varargin) == 1
    errorTh     = varargin{1};
elseif length(varargin) == 3
    errorTh     = varargin{1};
    lowerBound  = varargin{2};
    upperBound  = varargin{3};
elseif length(varargin) == 5
    errorTh     = varargin{1};
    lowerBound  = varargin{2};
    upperBound  = varargin{3};
    if strcmp(varargin{4},'Fit') && strcmp(varargin{5},'lsqcurvefit')
        lsqcurvefitFlag = 1;
    end
elseif length(varargin) == 7
    errorTh     = varargin{1};
    lowerBound  = varargin{2};
    upperBound  = varargin{3};
    if strcmp(varargin{4},'Fit') && strcmp(varargin{5},'lsqcurvefit')
        lsqcurvefitFlag = 1;
        lsqLowerBound = varargin{6};
        lsqUpperBound = varargin{7};
    end
elseif length(varargin) == 8
    errorTh     = varargin{1};
    lowerBound  = varargin{2};
    upperBound  = varargin{3};
    if strcmp(varargin{4},'Fit') && strcmp(varargin{5},'lsqcurvefit')
        lsqcurvefitFlag = 1;
        lsqLowerBound = varargin{6};
        lsqUpperBound = varargin{7};
    end
    display_flag = varargin{8};
else
    T1 = [];
    errordlg('Wrong list of arguments');
    return;
end
% END_initial -------------------------------------------------------------




%options for the fitting function
opts = optimset('Display', 'off');
%definition of the fit-model for mono-exponential fitting
mdl_exp_3 = @(beta,x)(  beta(1)*(1 - 2*exp(-x*beta(2)) + exp(-beta(2)*beta(3))  ) );


for i = 1:s1
    
    x = TIR(:)';
    y = bmCol(imagesTable(i, :))';
    N = size(x(:), 1); 
    
    
    % fitting -----------------------------------------------------
    
    
    z = max(y(:)) - y;
    z = z + max(z)/100; 
    
    [R1, A] = bmMonoExpFit(z, x, [], [], [], 'Fit', 'lsqcurvefit');  
    A       = A/2; 
    
    
    
    
    
    if isnan(A) || isinf(A) || isnan(R1) || isinf(R1)
        amp_map(i)                  = NaN;
        relax_map(i)                = NaN;
    else
        
        lsqLowerBound               = [0,       1/3000,     TR]; 
        lsqUpperBound               = [1000,    1/0.1,      TR]; 
        
        beta                        = [A, R1, TR];
        beta                        = lsqcurvefit(mdl_exp_3, beta, x, y, lsqLowerBound, lsqUpperBound, opts);
        amp_map(i)                  = beta(1);
        relax_map(i)                = beta(2);
        
        A                           = beta(1); 
        R1                          = beta(2); 
        
    end

    
if display_flag
        figure
        hold on
        plot(x, y, '.-')
        plot(x, A*(1 - 2*exp(-x*R1) + exp(-R1*TR)  ), 'x-');
end

end

if ndims(argImagesTable) > 2
    relax_map = reshape(relax_map, argSize(1, 1:end-1)); 
end

% END_fitting -------------------------------------------------


end

