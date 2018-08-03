function amat=funAlsopDP_vn(para)

nn=para.nn;
ranp=para.ranp;
IsAlsop=para.IsAlsop;
TD=para.TD;
df=para.df;

[Atr,Btr] = freeprecess(TD/2,para.T1,para.T2,df);

amat0=zeros(3,5);
amat0(3,1)=1;
amat0(:,2)=yrot(pi/2)*amat0(:,1); % the first 90 RF

amat0(:,3)=xrot(pi)*(Atr*amat0(:,2)+Btr); % relaxation and the 180 RF
amat0(:,4)=Atr*amat0(:,3)+Btr; % relaxation. before Alsop 
%amat0(:,5)=yrot(-pi/2)*amat0(:,4);

amat=repmat(amat0(:,4),1,nn);

% Apply Alsop gradient/Apply phase. 
if IsAlsop
    for ii=1:nn
        tmp=amat0(1,4).*exp((2*pi*ii/nn+ranp)*1i);
        amat(1,ii)=real(tmp);
        amat(2,ii)=imag(tmp);
    end
else
    tmp=amat0(1,4).*exp(ranp*1i);
    amat(1,:)=real(tmp);
    amat(2,:)=imag(tmp);
end

%Apply 90-RF & spoiler
for ii=1:nn
    amat(:,ii)=yrot(-pi/2)*amat(:,ii);
end
amat(1,:)=0; amat(2,:)=0;

end