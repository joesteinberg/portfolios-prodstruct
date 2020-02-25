function varargout = modeval(varargin);

nx=nargout;

na=length(varargin);

model=varargin{1};
eval(model{1});

if na>1
    i=2;
    while i<na
       pn=varargin{i};
       pv=varargin{i+1};
       tmp=strcat(pn,'=',num2str(pv),';');
       eval(tmp);
       i=i+2;
    end   
end

eval(model{2});
eval(model{3});

nd=model{4};
ns=model{5};
nb=model{6};
ne=model{7};
nu=model{8};

if nx==8  
    varargout{1}=A1;
    varargout{2}=A2;
    varargout{3}=A3;
    varargout{4}=A4;
    varargout{5}=B1;
    varargout{6}=B2;
    varargout{7}=SC;
    varargout{8}=SIGMA;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nx==7
    
    [dm nl]=size(A3);
    nr=nb+nu;
    nx=nd*(nd+1)/2;
    ny=ns*(ns+1)/2;

    AA1=[A1(1:nb,1:nb), A1(1:nb,ns+1:nd); A1(ns+1:nd,1:nb), A1(ns+1:nd,ns+1:nd)];
    AA2=[A2(1:nb,1:nb), A2(1:nb,ns+1:nd); A2(ns+1:nd,1:nb), A2(ns+1:nd,ns+1:nd)];

    AX1=[-A1(1:nb,nb+1:ns); -A1(ns+1:nd,nb+1:ns)];
    AX2=[A2(1:nb,nb+1:ns); A2(ns+1:nd,nb+1:ns)];

    NN=A2(nb+1:ns,nb+1:ns);

    AA3=AX1*NN+AX2;

    if SC(1,1)==-1
    
        AA4=0;
        AA5=0;
    
    else
    
        GG=zeros(nl,nx);
        GP=zeros(nl,nx);

        [nz, dm]=size(SC);
    
        for iz=1:nz
            i1=SC(iz,1);
            i2=SC(iz,2);
            i3=SC(iz,3);
            i4=SC(iz,4);
            cf=SC(iz,5);
            if i3<=ns;
                if i3<=nb;
                    i3=ne+i3;
                else
                    i3=i3-nb;
                end
            end
            if i4<=ns;
                if i4<=nb;
                    i4=ne+i4;
                else
                    i4=i4-nb;
                end
            end  
            if i3<=i4
                ic=(i4*(i4-1)/2)+i3;
            else
                tmp=i3;
                i3=i4;
                i4=tmp;
                ic=(i4*(i4-1)/2)+i3;
            end    
            if i2==0;
                GG(i1,ic)=cf;
            end    
            if i2==1;
                GP(i1,ic)=cf;
            end    
        end 

        AA4=[A3(1:nb,:); A3(ns+1:nd,:)]*GG;
        AA5=[A3(1:nb,:); A3(ns+1:nd,:)]*GP;

    end
    
    varargout{1}=AA1;
    varargout{2}=AA2;
    varargout{3}=AA3;
    varargout{4}=AA4;
    varargout{5}=AA5;
    varargout{6}=NN;
    varargout{7}=SIGMA;
    
end 
