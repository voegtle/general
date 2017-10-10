function [ N ] = hermgen( n, vector )
%hermgen HermGen generates a matrix which contains the values for n hermite
%polynomials at values designated by an input vector x
%   Using the recursion relation which defines Hermite polynomials, Hermgen
%   generates n polynomials and evaluates those polynomials at values
%   corresponding to those contained in the input vector x
% syms t 
nvect = (0:n);
[XX,NN] = meshgrid(vector,nvect);
%conditional case where n is either 0 or 1; program outputs row of ones
if n <= 1
    N = ones(1,length(vector));
else
    N = zeros(length(nvect),length(vector));
    %Case n == 1 and n == 2 are a bit unique for the recursion relation, so
    %I've added code to explicitly deal with these situations. 
    %The rest ofthe else statement 
    %deals with all other n values using the recursion
    %defined there (for reference, it is the probabilists definition of the
    %hermite polynomial. A for loop is used for some convenience, but
    %vectorized code does most of the work. 
    for i = 1:length(nvect)
        if i == 1
            N(i,:) = NN(i,:).*XX(i,:) + 1;
        elseif i == 2
            N(i,:) = N(i-1,:).*XX(i,:) - NN(i-1).*N(i-1,:);
        else
            N(i,:) = N(i-1,:).*XX(i,:) - NN(i-1).*N(i-2,:);
        end
    
    end
%This segment generates hermite polynomials symbolically, but it doesn't
%have much use for the assignment
%     hvect = sym(zeros(n,1));
%     hvect(1) = 1;
%     for i = 1:(length(hvect)-1)
%    
%         hvect(i+1) = t.*hvect(i) - diff(hvect(i),t);
%     end
%    
end

