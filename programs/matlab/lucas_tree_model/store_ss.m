% set all current exo/endo states and controls to steady state values
for i=1:NCOUNTRIES
    s=num2str(i);
    assignin('caller',['y_',s],0);
    assignin('caller',['c_',s],0);
    assignin('caller',['R_',s],log(1/beta));
    assignin('caller',['d_',s],log(theta));
    assignin('caller',['q_',s],log(beta*theta/(1-beta)));
    assignin('caller',['qm_',s],log(beta*theta/(1-beta)));
    assignin('caller',['WEALTH_',s],0);
end
for i=1:(NCOUNTRIES-1)
    assignin('caller',['Rx_',s,'3'],0);
end
for i=2:NCOUNTRIES
    assignin('caller',['cD_1',s],0);
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
