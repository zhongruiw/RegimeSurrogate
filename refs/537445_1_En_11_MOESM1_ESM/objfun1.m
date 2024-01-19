function f = objfun1(a)
global mean_truth variance_truth
% integration domain
integral_range = linspace(mean_truth - 4*sqrt(variance_truth), mean_truth + 4*sqrt(variance_truth), 200);
% loss function: entropy
f = trapz(integral_range, (a(1) + a(2) * integral_range + a(3) * integral_range.^2) .* exp(a(1) + a(2) * integral_range + a(3) * integral_range.^2));