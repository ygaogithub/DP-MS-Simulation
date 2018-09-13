function [MM, cmat]=funAlsopBssfp_vn(amat, para)
T1=para.T1;
T2=para.T2;
TRg=para.TRg;
TD=para.TD;
df=para.df;
RFthrad=para.RFthrad;
N_cata=para.N_cata;
necho=para.necho;
N_etl=para.N_etl;
nn=para.nn;
phi=para.phi;
IsAlsop=para.IsAlsop;
ranp=para.ranp;
if IsAlsop
    TR=para.TR+2*para.TAls;
    TE=para.TE+para.TAls;
else
    TR=para.TR;
    TE=para.TE;
end
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);
[Ate,Bte] = freeprecess(TE,T1,T2,df);

bmat=zeros(3,nn,necho);
cmat=zeros(3,nn,N_etl);
% start the bssfp simulation
for ii=1:nn
    n_count=1;
    if IsAlsop
        sig_rf=zrot(phi)*xrot(0.5*RFthrad/N_cata)*amat(:,ii);
        sig_alsop=funCalAlsopDephasing(sig_rf,-2*pi*ii/nn);
        sig_new=Atr*(Ate*sig_alsop+Bte)+Btr;
        sig_alsop2=funCalAlsopDephasing(sig_new,2*pi*ii/nn);
        bmat(:,ii,1)=sig_alsop2;
    else
        bmat(:,ii,1) = Atr*(Ate*zrot(phi)*xrot(0.5*RFthrad/N_cata)*amat(:,ii)+Bte)+Btr;
    end
    
    % The following RF pulse of the cata
    for nprep=2:N_cata
        Rflip_prep=zrot(phi)*xrot((nprep-0.5)*RFthrad/N_cata);
        if IsAlsop
            sig_rf=Rflip_prep*bmat(:,ii,nprep-1);
            sig_alsop=funCalAlsopDephasing(sig_rf,-2*pi*ii/nn);
            sig_new=Atr*(Ate*sig_alsop+Bte)+Btr;
            sig_alsop2=funCalAlsopDephasing(sig_new,2*pi*ii/nn);
            bmat(:,ii,nprep)=sig_alsop2;
        else
            bmat(:,ii,nprep)=Atr*(Ate*Rflip_prep*bmat(:,ii,nprep-1)+Bte)+Btr;
        end
    end
end

for ii=1:nn
    % Pulses of bssfp seq
    Rflip = zrot(phi)*xrot(RFthrad);
    for n_ro=1:N_etl
        n_count=n_count+1;
        sig_rf=Rflip*bmat(:,ii,n_ro+nprep-1);
        if IsAlsop
            sig_alsop=funCalAlsopDephasing(sig_rf,-2*pi*ii/nn);
            sig_TE=Ate*sig_alsop+Bte;
        else
            sig_TE=Ate*sig_rf+Bte;
        end
        cmat(:,ii,n_ro)=sig_TE;
        
        if IsAlsop
            sig_TR=Atr*sig_TE+Btr;
            sig_alsop2=funCalAlsopDephasing(sig_TR,2*pi*ii/nn);
            bmat(:,ii,n_ro+nprep)=sig_alsop2;
        else
            bmat(:,ii,n_ro+nprep)=Atr*sig_TE+Btr;
        end
        
    end
end

MM=squeeze(sum(cmat,2));
% Mxy= sum((squeeze(cmat(1,:,:)+1i*cmat(2,:,:))),1);
% Mxy_ampl= abs(Mxy);
% Mxy_angle=angle(Mxy);

%relaxization 
TR_relax=TRg-TD-TR*(N_etl+N_cata);
[Afpl,Bfpl] = freeprecess(TR_relax,T1,T2,df);
MM=MM./nn;
MM(:,N_etl+1) = Afpl*MM(:,N_etl)+Bfpl;
Mxy= sum((squeeze(MM(1,:,:)+1i*MM(2,:,:))),1);
Mxy_ampl= abs(Mxy);
Mxy_angle=angle(Mxy);

if para.fplot
    if IsAlsop
        fname='With Alsop Dephasing';
    else
        fname='Without Alsop Dephasing';
    end
    
    figure;
    subplot(1,2,1); plot(Mxy_ampl,'-*');xlabel('TR #');
    ylabel('Mxy magnitude New');  title([fname, ', phase= ' num2str(ranp/pi)]);
    
    hold on;
    subplot(1,2,2); plot(Mxy_angle,'-*');xlabel('TR #');
    ylabel('Mxy phase');title([fname, ', df= ' num2str(df)]);
end

end