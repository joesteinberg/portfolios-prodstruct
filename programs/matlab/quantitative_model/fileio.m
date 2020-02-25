inin_base = zeros(NCOUNTRIES,NCOUNTRIES);
va_base = zeros(NCOUNTRIES,3);
y_base = zeros(NCOUNTRIES,1);
ff_base = zeros(NCOUNTRIES,2,NCOUNTRIES);

inin_raw = importdata([input_path, 'inin',str,mname,'.txt']);
va_raw = importdata([input_path, 'va',str,mname,'.txt']);
fin_raw = importdata([input_path, 'fin',str,mname,'.txt']);

for row=1:size(inin_raw.data,1)
    i=inin_raw.data(row,1);
    %s=inin_raw.data(row,2);
    j=inin_raw.data(row,3);
    %r=inin_raw.data(row,4);
    val=inin_raw.data(row,5);
    inin_base(i,j)=val;        
end

for row=1:size(va_raw.data,1)
    i=va_raw.data(row,1);
    %s=va_raw.data(row,2);
    val1=va_raw.data(row,3);
    val2=va_raw.data(row,4);
    y_base(i) = val1;
    va_base(i,3) = val2;
end

for row=1:size(fin_raw.data,1)
    i=fin_raw.data(row,1);
    j=fin_raw.data(row,2);
    %s=fin_raw.data(row,3);
    val1=fin_raw.data(row,4);
    val2=fin_raw.data(row,5);
    ff_base(i,1,j) = val1;
    ff_base(i,2,j) = val2;
end
