%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% doing one-time setup stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str='-';

fprintf('\n//////////////////////////////////////////////////////////////////\n')
fprintf(['INTERNATIONAL PORTFOLIO DIVERSIFICATION AND THE STRUCTURE ' ...
         'OF GLOBAL PRODUCTION\nJoseph Steinberg, University of ' ...
         'Toronto\nLast Updated: June, 2017\n']);

create_model=1;
if create_model == 1

    fprintf(['\n/////////////////////////////////////////////////////' ...
             '/////////////\n'])
    fprintf('PERFORMING ONE-TIME SETUP...\n\n');

    NCOUNTRIES=4;
    TINY=1.0e-9;
    Tsim=10010;

    addpath '../usg_toolkit';
    addpath '../devereux_sutherland_example';

    input_path = '../../python/calibration_data/';
    output_path = 'output/';
    cnames = {'USA','ADV','EME','ROW'};
    mnames = {'bench-1995',...
              'bench-2011',...
              'size-counter',...
              'trd-counter',...
              'io-counter',...
              'nx-counter',...
              'bal-bench-1995',...
              'bal-bench-2011',...
              'bal-size-counter',...
              'bal-trd-counter',...
              'bal-io-counter'};

    mtitles = {'1995 benchmark',...
               '2011 benchmak',...
               'Changes in size only',...
               'Changes in openness only',...
               'Changes in intermediate trade only',...
               'Changes in trade balances only',...
               'Balanced-trade 1995 benchmark',...
               'Balanced-trade 2011 benchmark',...
               'Balanced-trade size counterfactual'...
               'Balanced-trade openness counterfactual',...
               'Balanced-trade IO counterfactual'};

    fprintf(['\tDeclaring symbolic variables and constructing analytical ' ...
             'representation of equilibrium conditions\n']);
    tic;
    analytical_params_vars_eqns;
    fprintf('\t');
    toc;

    fprintf(['\n\tTaking analytical derivatives of equilibrium system\n']);
    tic;
    analytical_derivs;
    fprintf('\t');
    toc;
end

fprintf(['\n/////////////////////////////////////////////////////' ...
         '/////////////\n'])
%tic;
%fprintf(['Loading stochastic process parameters and drawing shocks for simulation\n']);
%SIGMA=importdata([input_path,'SIGMA.txt']);
%N=importdata([input_path,'N.txt']);
%EPSsim = mvnrnd(zeros(1,nz),SIGMA,Tsim)';
%toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calibrate and solve model using benchmark elasticities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose=0;
sens=1;
noio=0;

% vector of wedges... we choose this when
% solving the 1995 model so that we get the portfolio
% diversification figures exactly right in that year
% also store the targets for this part of the calibration in DIV
NT=4;
TAU = zeros(NT,1);
TAU0 = TAU;
TAU1 = TAU;
da = importdata('../../python/results/div_for_matlab.csv');
divcol=4;
DIV = da.data(1:NT,divcol);

for sens=0:(8*sens)

    fprintf('\n//////////////////////////////////////////////////////////////////\n')
   
    if(sens==0)
        fprintf('BASELINE MODEL\n');
    elseif(sens==1)
        fprintf('HIGHER ARMINGTON ELASTICITY\n')
    elseif(sens==2)
        fprintf('LOWER ARMINGTON ELASTICITY\n')
    elseif(sens==3)
        fprintf('HIGHER FINAL ARMINGTON ELASTICITY\n');
    elseif(sens==4)
        fprintf('LOWER FINAL ARMINGTON ELASTICITY\n');
    elseif(sens==5)
        fprintf('COBB-DOUGLAS PRODUCTION\n');
    elseif(sens==6)
        fprintf('UNCORRELATED SHOCKS\n');
    elseif(sens==7)
        fprintf('UNCORRELATED SHOCKS + COBB-DOUGLAS PRODUCTION\n');
    elseif(sens==8)
        fprintf('HIGHER RISK AVERSION\n');
    end

    if(sens==6||sens==7)
        %tic;
        %fprintf(['\tReloading stochastic process parameters and redrawing shocks for simulation\n']);
        SIGMA=importdata([input_path,'SIG2.txt']);
        N=importdata([input_path,'N2.txt']);
        %EPSsim = mvnrnd(zeros(1,nz),SIGMA,Tsim)';
        %toc;
    else
        SIGMA=importdata([input_path,'SIGMA.txt']);
        N=importdata([input_path,'N.txt']);
    end

    for outer=1:size(mnames,2)
        tic

        if(outer==1 || outer==7)
            TAU=zeros(NT,1);
            tauflag=1;
        else
            tauflag=0;
        end

        mname = mnames{outer};
        mtitle = mtitles{outer};

        fprintf('\n\t------------------------------------------\n')
        fprintf(['\tScenario: ',mnames{outer},'\n\n']);

        fprintf(['\t\tReading input-output matrix data\n']);
        fileio;

        fprintf(['\t\tCalibrating model parameters\n']);
        calibration;

        if(sens==0 & outer==7)
            %calibration_latex1;
        end

        fprintf(['\t\tStoring numeric steady state in symbolic variables\n']);
        store_ss;

        fprintf(['\t\tEvaluating approximated equilibrium conditions numerically\n']);
        approx=1;
        for j=1:length(FF)
            tmp = eval(FF(j));
            if(abs(tmp)>1.0e-5)
                error(['\t\tEquilibrium condition ',num2str(j),...
                       ' does not hold in steady state! Value = ',num2str(tmp)]);
            end
        end

        fprintf(['\t\tSolving for steady state portfolios\n']);
        if(tauflag==1)
            fprintf(['\t\t1995 benchmark... calibrating wedges\n']);
        end
        zero_order_portfolios;

        %if sens==0
        %    fprintf(['\t\tSimulation\n']);
        %    simul;
        %end
        
        fprintf('\t\t');
        toc 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create latex tables for main results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n//////////////////////////////////////////////////////////////////\n')
fprintf(['WRITING RESULTS TO LATEX FILES\n']);
calibration_latex2;
results_latex;
