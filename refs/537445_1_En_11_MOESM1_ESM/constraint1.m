function [c,ceq]=constraint1(a)
global mean_truth variance_truth
% integration interval
integral_range = linspace(mean_truth - 4*sqrt(variance_truth), mean_truth + 4*sqrt(variance_truth), 200);

% the three constaints
c0 = trapz(integral_range, exp(a(1) + a(2) * integral_range + a(3) * integral_range.^2)) - 1;
c1 = trapz(integral_range, integral_range .* exp(a(1) + a(2) * integral_range + a(3) * integral_range.^2)) - mean_truth;
c2 = trapz(integral_range, (integral_range-c1).^2 .* exp(a(1) + a(2) * integral_range + a(3) * integral_range.^2)) - variance_truth;
% equality and inequality constraints; the equality constraints are
% represented by the inequalities for the purpose of numerical stability
c = [];
ceq = [c0;c1;c2;-c0;-c1;-c2];