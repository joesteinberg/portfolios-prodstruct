function [F,ALPHA,ALPHA3,lambda] = zero_order_portfolios2(TAU,params)
    NT=4;
    beta = params.beta;
    WEALTH_1 = params.WEALTH_1;
    WEALTH_2 = params.WEALTH_2;
    WEALTH_3 = params.WEALTH_3;
    WEALTH_4 = params.WEALTH_4;
    q_1 = params.q_1;
    q_2 = params.q_2;
    q_3 = params.q_3;
    q_4 = params.q_4;
    NCOUNTRIES = params.NCOUNTRIES;
    ALPHA0 = params.ALPHA0;
    TINY = params.TINY;
    DIV = params.DIV;
    SIGMA = params.SIGMA;

    % ----------------------------------------------------------------------
    % solve for steady-state portfolios
    % ----------------------------------------------------------------------

    options = optimset('fsolve');
    options = optimset(options,'TolFun',1e-8);
    options = optimset(options,'TolX',0.0);
    options = optimset(options,'MaxIter',10000);
    options = optimset(options,'MaxFunEvals',40000);
    options = optimset(options,'Display','off');
    options = optimset(options,'Algorithm','trust-region-dogleg');
    
    TAU2 =TAU;
    [ALPHA,fval,exitflag] = fsolve( @(ALPHA) share_func(ALPHA,params,TAU2) , ALPHA0, options);
    if(exitflag ~= 1 && max(abs(fval/1000)>1e-12))
        error('Warning! Fsolve failed to find equilibrium portfolio shares!');
    end

    % calculate full matrix of asset holdings
    ALPHA2 = [ALPHA*beta,zeros(NCOUNTRIES-1,1);zeros(1,NCOUNTRIES)];
    for i=1:(NCOUNTRIES-1)
        s=num2str(i);
        ALPHA2(i,NCOUNTRIES) = eval(['WEALTH_',s]) - sum(ALPHA2(i,:));
    end
    ALPHA2(NCOUNTRIES,:) = -sum(ALPHA2(1:(NCOUNTRIES-1),:),1);

    % check that everything adds up
    for i=1:NCOUNTRIES
        s=num2str(i);
        tmp = eval(['WEALTH_',s]) - sum(ALPHA2(i,:));
        if abs(tmp)>TINY
            error(['alpha2(',s,':) do not sum to WEALTH(',s,')!']);
        end
        tmp = sum(ALPHA2(:,i));
        if abs(tmp)>TINY
            error(['alpha2(',s,':) do not sum to zero!']);
        end
    end

    % calculate portfolio weights
    ALPHA3 = ALPHA2;
    for i=1:NCOUNTRIES
        s=num2str(i);
        ALPHA3(i,i) = ALPHA3(i,i) + exp(eval(['q_',s]));
        ALPHA3(i,:) = ALPHA3(i,:)/(eval(['WEALTH_',s]) + exp(eval(['q_',s])));
    end

    % calculate lambda (holdings of actual shares)
    lambda = zeros(size(ALPHA2));
    for i=1:(NCOUNTRIES)
        for j=1:NCOUNTRIES
            tmp = ALPHA2(i,j);
            if(i==j)
                tmp = tmp + exp(eval(['q_',num2str(i)]));
            end
            lambda(i,j) = tmp/exp(eval(['q_',num2str(j)]));
        end
    end

    % calculate errors
    F = DIV(:) - 100*(1-diag(ALPHA3(1:NT,1:NT)));
end
