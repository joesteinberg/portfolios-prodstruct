for i=1:NCOUNTRIES
    s=num2str(i);

    assignin('caller',['WEALTH_',s],b_base(i));

    assignin('caller',['ngdp_',s],log(GDP_base(i)))
    assignin('caller',['rgdp_',s],log(GDP_base(i)))
    assignin('caller',['y_',s],log(y_base(i,1)))
    assignin('caller',['va_',s],log(va_base(i,3)))
    assignin('caller',['m_',s],log(M_base(i)))
    assignin('caller',['fc_',s],log(f_base(i,1)))
    assignin('caller',['fi_',s],log(f_base(i,2)))
    assignin('caller',['k_',s],log(K_base(i)))
    assignin('caller',['n_',s],log(L_base(i)))
    assignin('caller',['ex_',s],log(EX_base(i)))
    assignin('caller',['IM_',s],log(IM_base(i)))
    assignin('caller',['NX_',s],log(NX_base(i)))
    assignin('caller',['div_',s],y_base(i) - L_base(i) - f_base(i,2) ...
             - M_base(i));

    for j=1:NCOUNTRIES
        s2 = num2str(j);
        assignin('caller',['m_',s,s2],log(inin_base(i,j)));
        assignin('caller',['fc_',s,s2],log(ff_base(i,1,j)));
        assignin('caller',['fi_',s,s2],log(ff_base(i,2,j)));
    end

    assignin('caller',['R_',s],log(1/beta));
    assignin('caller',['rk_',s],log(r_base+ eval(['delta_',s])));
    assignin('caller',['q_',s],log(beta*eval(['div_',s])/(1-beta)));
    assignin('caller',['qm_',s],log(beta*eval(['div_',s])/(1-beta)));
    assignin('caller',['pfc_',s],0.0);
    assignin('caller',['pfi_',s],0.0);
    assignin('caller',['py_',s],0.0);
    assignin('caller',['w_',s],0.0);
    assignin('caller',['tfp_',s],0.0);

    if i>1
        assignin('caller',['rer_1',s],0.0);
    end
    if i<NCOUNTRIES
        assignin('caller',['cD_',s,num2str(NCOUNTRIES)],0.0);
        assignin('caller',['Rx_',s,num2str(NCOUNTRIES)],0.0);
    end

    assignin('caller',['wn_',s],log(L_base(i)));
    assignin('caller',['qd_',s],eval(['div_',s])+ ...
                                    exp(eval(['q_',s])));
end

% next period
for i=1:length(ZZ)
    name = char(ZZ(i));
    assignin('caller',[name,'_p'],eval(name));
end

for i=1:length(XX)
    name = char(XX(i));
    assignin('caller',[name,'_p'],eval(name));
end

for i=1:length(YY)
    name = char(YY(i));
    assignin('caller',[name,'_p'],eval(name));
end

