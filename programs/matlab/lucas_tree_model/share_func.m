function F = share_func(alphatt,params)
    NCOUNTRIES = params.NCOUNTRIES;
    R1 = params.R1;
    R2 = params.R2;
    D1 = params.D1;
    D2 = params.D2;
    SIGMA = params.SIGMA;

    H = inv(eye(NCOUNTRIES-1)-alphatt*R1)*(alphatt*R2);
    Rt = R1*H+R2;
    Rts = Rt*SIGMA;

    F = [];
    for i=1:(NCOUNTRIES-1)
        Dt = D1(i,:)*H + D2(i,:);
        F=[F;Rts*(Dt')];
    end
    F = 1000*F;
end