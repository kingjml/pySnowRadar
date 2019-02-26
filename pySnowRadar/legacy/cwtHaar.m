function coefs = cwtHaar(x,scales)
% cwtHaar - returns the continuous wavelet transform (CWT) of the 
% real-valued signal S. The wavelet transform is computed for the specified 
% odd number scales using the analyzing wavelet Haar. 
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    x      - The input data
%    scales - A 1-D vector with positive elements. All scales must be odd
%             and monotonically increasing.
%
% Outputs:
%    coefs - A matrix with the number of rows equal to the length of scales 
%            and number of columns equal to the length of the input signal. 
%            The k-th row of coefs corresponds to the CWT coefficients for 
%            the k-th elemens in the scales vector.
%
% Example:
% -
%
% Other m-files required: none
% 
% Subfunctions: none
%
% MAT-files required: none
%
% See also: waveletSnowProcessor
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
%
% August 2016; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if any(mod(scales,2)==0)
   error('Error. Scales must all be odd')
end
% Check to see if there are any even scales.

dataRow = x(:).';
% Make row vector

maxHalfScale = floor(scales(end)/2);
% The maximum half scale.


upperTri = bsxfun(@rdivide,triu(ones(maxHalfScale),0),sqrt(scales));
% The upper triangle.

zeroVect = zeros(1,maxHalfScale);
% The middle zero vector.

triComb = [-flipud(upperTri);zeroVect;upperTri];
% Combining the triangular elements.

triCombT = triComb.';
% The transpose of the tricomb.


C = conv2(1,dataRow,triCombT,'full');
% Perform 2D convolution.

startCol = maxHalfScale+1;
endCol   = length(dataRow)+maxHalfScale;
% The column limits containing the data.

coefs = C(:,startCol:endCol);
% The wavelet coeficients. 

%==========================================================================
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<