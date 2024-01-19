% Chapter 2
% Given the first two or four moments, finding the least biased PDF based
% on the maximum entropy principle.
% Input: Moment_input, which is a 2- or 4-dimensional vector containing the
% first two or the first four moments, namely [mean, variance] or [mean,
% variance, skewness, kurtosis].
% Note that the skewness and kurtosis cannot be arbitary since otherwise
% the associated PDF may not exist
% Example of input: 
% Moment_input = [0,1];
% Moment_input = [0,1,1,4];
% Output: the coefficients 'Coeff' of the exponential function  
% exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2), or 
% exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4)
% which is the parametric form of the PDF
function Coeff = Max_Entropy_Principle_withMoments(Moment_input)
num_of_moments = length(Moment_input); % detect the length of the input 
% define the global variables of the truth of the moments 
global mean_truth variance_truth skewness_truth kurtosis_truth
if num_of_moments == 2 % if the input only contains mean and variance
    mean_truth = Moment_input(1); variance_truth = Moment_input(2);
    [Coeff,~]=fmincon('objfun1',[0,0,0],[],[],[],[],[],[],'constraint1');
elseif num_of_moments == 4 % if the input contains the leading four moments
    mean_truth = Moment_input(1); variance_truth = Moment_input(2);
    skewness_truth = Moment_input(3); kurtosis_truth = Moment_input(4);
    [Coeff,~]=fmincon('objfun2',[0,0,0,0,0],[],[],[],[],[],[],'constraint2');
else % if the input has a different length than 2 or 4, then quit the code
    disp('Error')
    disp('Need to provide [mean, variance] or [mean, variance, skewness, kurtosis]')
    return
end
    
% integration interval: mean plus/minus four standard deviation 
integral_range = linspace(mean_truth - 4*sqrt(variance_truth), mean_truth + 4*sqrt(variance_truth), 200);


% showing the results of the moments computed from the reconstructed PDF and comparing with the truth 
if num_of_moments == 2
    normalization_estimation = trapz(integral_range, exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2));
    mean_estimation          = trapz(integral_range, integral_range .* exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2));
    variance_estimation      = trapz(integral_range, (integral_range - mean_estimation).^2 .* exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2));
    disp('The true value:')
    disp('normalization constant; mean; variance')
    disp([1, mean_truth, variance_truth])
    disp('Using the estimated PDF from maximum entropy principle:')
    disp('normalization constant; mean; variance')
    disp([normalization_estimation, mean_estimation, variance_estimation])
end
if num_of_moments == 4
    normalization_estimation = trapz(integral_range, exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4));
    mean_estimation          = trapz(integral_range, integral_range .* exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4));
    variance_estimation      = trapz(integral_range, (integral_range - mean_estimation).^2 .* exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4));
    skewness_estimation      = trapz(integral_range, (integral_range - mean_estimation).^3 .* exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4))/variance_estimation^(3/2);
    kurtosis_estimation      = trapz(integral_range, (integral_range - mean_estimation).^4 .* exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4))/variance_estimation^2;
    disp('The true value:')
    disp('normalization constant; mean; variance, skewness, kurtosis')
    disp([1, mean_truth, variance_truth, skewness_truth, kurtosis_truth])
    disp('Using the estimated PDF from maximum entropy principle:')
    disp('normalization constant; mean; variance, skewness, kurtosis')
    disp([normalization_estimation, mean_estimation, variance_estimation, skewness_estimation, kurtosis_estimation])
end
% showing the reconstructed PDF
figure
if num_of_moments == 2
    plot(integral_range, exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2),'b','linewidth',2)
elseif num_of_moments == 4
    plot(integral_range, exp(Coeff(1) + Coeff(2) * integral_range + Coeff(3) * integral_range.^2 + Coeff(4) * integral_range.^3 + Coeff(5) * integral_range.^4),'b','linewidth',2)
end
box on
set(gca,'fontsize',12)
title('Reconstructed PDF via maximum entropy principle')
xlabel('x')
ylabel('Probability')