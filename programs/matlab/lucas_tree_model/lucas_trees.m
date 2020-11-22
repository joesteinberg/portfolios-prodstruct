fprintf(['\nINTERNATIONAL PRODUCTION CHAINS AND THE PORTFOLIO ' ...
         'DIVERSIFICATION PUZZLE\n']);
fprintf(['Three-country Lucas tree example model\n\n']);

NCOUNTRIES=7;
TINY=1.0e-9;
addpath '../usg_toolkit';
addpath '../devereux_sutherland_example';

fprintf(['Declaring symbolic variables and constructing analytical ' ...
         'representation of equilibrium conditions\n']);
tic;
analytical;
toc;

fprintf(['\nTaking analytical derivatives of equilibrium system\n']);
tic;
analytical_derivs;
toc;

gama = 1;
theta = 1/3;
beta = 1/1.04;
N = 0.91*eye(NCOUNTRIES);
SIGMA = 0.06*eye(NCOUNTRIES);

fprintf(['\nStoring numeric steady state in symbolic variables\n']);
store_ss;

fprintf(['\nEvaluating approximated equilibrium conditions numerically\n']);
approx=1;
for i=1:length(FF)
    tmp = eval(FF(i));
    if(abs(tmp)>1.0e-5)
        error(['Equilibrium condition ',num2str(i),...
               ' does not hold in steady state! Value = ',num2str(tmp)]);
    end
end

fprintf(['\nSolving for steady state portfolios\n']);
zero_order_portfolios;
