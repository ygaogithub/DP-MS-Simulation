%Simulate nn spins with nn different phase
nn=2000;
para.nn=nn;
para.TR=3.7;
para.TAls=0.5; % alsop gradient duration 
para.TD=40;
para.TE=para.TR/2;
para.RFthrad=pi/3; % flip angle
para.N_cata=10;
para.phi=pi;
para.T1 = 800;
para.T2 = 80;
para.df=0;
para.N_etl=48;
para.necho=para.N_cata+para.N_etl;
para.IsAlsop=1;
para.fplot=0;
para.ranp=0; % random added phase

%% Mag/phase with/without Alsop
np=31; % number of simulated phase. 
Mxy1=zeros(np,para.N_etl);
Mxy2=Mxy1;
rphase=linspace(-pi,pi,np);
for ii=1:np
    % without Alsop
    para.ranp=rphase(ii);
    para.IsAlsop=0;
    amat=funAlsopDP_vn(para);
    Mxy1(ii,:)=funAlsopBssfp_vn(amat, para);
    
    % with Alsop
    para.ranp=rphase(ii);
    para.IsAlsop=1;
    amat=funAlsopDP_vn(para);
    Mxy2(ii,:)=funAlsopBssfp_vn(amat, para);
end
%%
nline=1;
figure;
plot(rphase,abs(Mxy1(:,nline)),'b-o', rphase,abs(Mxy2(:,nline)),'r-*'); title('Mag')


