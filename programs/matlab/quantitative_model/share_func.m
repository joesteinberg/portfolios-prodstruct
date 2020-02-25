function F = share_func(alphatt,params,TAU2)
    NCOUNTRIES = params.NCOUNTRIES;
    R1 = params.R1;
    R2 = params.R2;
    D1 = params.D1;
    D2 = params.D2;
    SIGMA = params.SIGMA;
    WEALTH_1 = params.WEALTH_1;
    WEALTH_2 = params.WEALTH_2;
    WEALTH_3 = params.WEALTH_3;
    WEALTH_4 = params.WEALTH_4;
    beta = params.beta;

    %a2 = [alphatt*beta,zeros(NCOUNTRIES-1,1);zeros(1,NCOUNTRIES)];
    %for i=1:(NCOUNTRIES-1)
    %    s=num2str(i);
    %    a2(i,NCOUNTRIES) = eval(['WEALTH_',s]) - sum(a2(i,:));
    %end
    %a2(NCOUNTRIES,:) = -sum(a2(1:(NCOUNTRIES-1),:),1);

    H = inv(eye(NCOUNTRIES-1)-alphatt*R1)*(alphatt*R2);
    Rt = R1*H+R2;
    Rts = Rt*SIGMA;

    F = [];
    for i=1:(NCOUNTRIES-1)
        Dt = D1(i,:)*H + D2(i,:);
        tmp = (Rts)*(Dt');
        for j=1:(NCOUNTRIES-1)
            tmp(j) = tmp(j) + 1.0 - exp(TAU2(NCOUNTRIES)) - exp(TAU2(i));
            %tmp(j) = tmp(j) + TAU(NCOUNTRIES)*(a2(NCOUNTRIES,j)^2.0) + TAU(i)*(a2(i,NCOUNTRIES)^2.0);
            if(i ~=j)
                tmp(j) = tmp(j) + exp(TAU2(i));
                %tmp(j) = tmp(j) - TAU(i)*(a2(i,j)^2.0);
            else
                tmp(j) = tmp(j) + 1.0;
                %tmp(j) = tmp(j) + 0.0;
            end
        end
        F=[F;tmp];
    end

    F=1000*F;
end