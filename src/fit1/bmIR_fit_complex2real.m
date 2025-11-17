% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023


function b = bmIR_fit_complex2real(a)

argSize = size(a); 
s1      = prod(bmCol(  argSize(1, 1:end-1)  )); 
s2      = argSize(1, end);  

a = reshape(a, [s1, s2]); 

phase_1     = angle(a(:, 1));
phase_1     = repmat(phase_1,   [1, s2]); 

phase_end   = angle(a(:, end));
phase_end   = repmat(phase_end, [1, s2]); 

diff_1      = abs(angle(a) - phase_1); 
diff_end    = abs(angle(a) - phase_end); 


myMin       = min(diff_1, diff_end); 

m_1         = (myMin == diff_1);
m_end       = (myMin == diff_end);

b           = abs(a).*m_end - abs(a).*m_1; 

[~, minInd] = min(abs(b), [], 2); 
minInd      = repmat(minInd, [1, s2]); 

myInd       = 1:s2; 
myInd       = repmat(myInd, [s1, 1]); 

m_neg       = ((myInd - minInd) < 0);
m_0         = ((myInd - minInd) == 0);
m_pos       = ((myInd - minInd) > 0); 


b           = abs(b).*m_pos + b.*m_0 - abs(b).*m_neg;  

b           = reshape(b, argSize); 

end