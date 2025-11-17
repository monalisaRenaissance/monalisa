% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023


function [a_map b_map c_map varargout] = bmComplexMonoExpFit1(argImagesTable, argX, varargin)

mySize = size(argImagesTable);
mySize = [prod(mySize(1:end-1)) mySize(end)];

if not(length(argX) == mySize(2))
    a_map = 0;
    b_map = 0;
    c_map = 0;
    d_map = 0;
    errordlg('Wrong list of arguments');
    return;
end

errorTh = [];
lowerBound = [];
upperBound = [];
lsqLowerBound = [];
lsqUpperBound = [];
lsqcurvefitFlag = 0;

if length(varargin) == 0
    1+1;
elseif length(varargin) == 1
    errorTh = varargin{1};
elseif length(varargin) == 3
    errorTh = varargin{1};
    lowerBound = varargin{2};
    upperBound = varargin{3};
else
    a_map = 0;
    b_map = 0;
    c_map = 0;
    d_map = 0;
    errordlg('Wrong list of arguments');
    return;
end


%definition of the fit-model for mono-exponential fitting
mdl_mono_exp = @(beta,x)(beta(1)*exp(-x*beta(2)));


%options for the fitting function
opts = optimoptions(@lsqcurvefit, 'Display', 'off', 'Algorithm', 'levenberg-marquardt', 'FunValCheck', 'off');

imagesTable = reshape(argImagesTable, mySize);
iMax = mySize(2);
x = squeeze(argX)';
x = reshape(x, [1 length(x)]);
xTable = repmat(x, [mySize(1) 1]);

a_map   = zeros(mySize(1), 1);
b_map   = zeros(mySize(1), 1);
c_map   = zeros(mySize(1), 1);
d_map   = zeros(mySize(1), 1);

rhoTable = abs(imagesTable);
alphaTable = angle(imagesTable);

[d_map, c_map, alphaFit] = bmAffinePhaseFit(alphaTable, x);
imagesTableShifted = imagesTable.*exp(-i*alphaFit);

logRhoTable = log(abs(imagesTable));

MeanX = mean(xTable, 2);
MeanZ = mean(logRhoTable, 2);
MeanX2 = mean(xTable.^2, 2);
MeanXZ = mean(xTable.*logRhoTable, 2);

h = (MeanX2.*MeanZ-MeanX.*MeanXZ)./(MeanX2-MeanX.^2);
a_map_start = exp(h);
b_map_start = -(MeanXZ-MeanX.*MeanZ)./(MeanX2-MeanX.^2);

a_map = a_map_start;
b_map = b_map_start;

for k = 1:mySize(1)
    if isnan(a_map_start(k))||isnan(b_map_start(k))
        a_map(k) = NaN;
        b_map(k) = NaN;
        c_map(k) = NaN;
    else
        y = imagesTableShifted(k, :);
        
        beta = [a_map_start(k) b_map_start(k)];
        beta = lsqcurvefit(mdl_mono_exp , beta, x, y, [], [], opts);
        a_map(k) = beta(1);
        b_map(k) = beta(2);
    end
end

a_map_table = repmat(a_map, [1 length(x)]);
b_map_table = repmat(b_map, [1 length(x)]);
c_map_table = repmat(c_map, [1 length(x)]);
d_map_table = repmat(d_map, [1 length(x)]);

myFit = a_map_table.*exp(-b_map_table.*xTable);
myFit = myFit.*exp(i*alphaFit);  
myError = sqrt(mean(abs(myFit-imagesTable).^2./abs(myFit).^2,2));

if not(isempty(errorTh))
    errorMask = (myError > errorTh);
else
    errorMask = zeros(mySize(1), 1);
end
errorMask = errorMask + isnan(a_map)+isnan(b_map)+isnan(c_map);
errorMask = logical(errorMask);

if not(isempty(lowerBound))
    errorMask = errorMask + (b_map < lowerBound);
    errorMask = logical(errorMask);
end
if not(isempty(upperBound))
    errorMask = errorMask + (b_map > upperBound);
    errorMask = logical(errorMask);
end

a_map(errorMask) = NaN;
b_map(errorMask) = NaN;
c_map(errorMask) = NaN;
d_map(errorMask) = NaN;


argSize = size(argImagesTable);
mySize = argSize(1:end-1);

if ndims(argImagesTable) > 2
    errorMask = reshape(errorMask, mySize);
    a_map = reshape(a_map, mySize);
    b_map = reshape(b_map, mySize);
    c_map = reshape(c_map, mySize);
end

% a_map = real(a_map);
% b_map = real(b_map); 

varargout{1} = reshape(myFit, argSize);
varargout{2} = errorMask;

end
