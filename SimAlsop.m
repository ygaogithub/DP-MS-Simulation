%Simulate nn spins with nn different phase
nn=2000;
para.nn=nn;
para.TR=3.7;
para.TRg=2000;
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
para.N_shot=5; % shot number. The first shot should be discarded due to signal not reaching steady level.
para.necho=para.N_cata+para.N_etl;
para.IsAlsop=1;
para.fplot=0;
para.ranp=0; % random added phase

%% Mag/phase with/without Alsop
np=31; % number of simulated phase.
MM1=zeros(3,para.N_etl+1,np,para.N_shot);
MM2=MM1;
rphase=linspace(-pi,pi,np);
for jj=1:para.N_shot
    
    for ii=1:np
        % without Alsop
        para.ranp=rphase(ii);
        para.IsAlsop=0;
        if jj==1;  amat0=[0 0 1];
        else  amat0=MM1(:,end,ii, jj-1);
        end
        amat=funAlsopDP_vn(amat0,para);
        MM1(:,:,ii,jj)=funAlsopBssfp_vn(amat, para);
        
        % with Alsop
        para.ranp=rphase(ii);
        para.IsAlsop=1;
        if jj==1;  amat0=[0 0 1];
        else  amat0=MM2(:,end, ii, jj-1);
        end
        amat=funAlsopDP_vn(amat0,para);
        MM2(:,:,ii,jj)=funAlsopBssfp_vn(amat, para);
    end
end
%%
necho=1;nshot=2;
figure;
Mxy1= squeeze(MM1(1,:,:,:)+1i*MM1(2,:,:,:));
Mxy2= squeeze(MM2(1,:,:,:)+1i*MM2(2,:,:,:));
plot(rphase,abs(Mxy1(necho,:,nshot)),'b-o', rphase,abs(Mxy2(necho,:,nshot)),'r-*');
set(gca,'FontSize',18); axis([-3.5 3.5 0 0.35]);
xlabel('Phase Error (radians)'); ylabel('a.u.'); legend('w/o MS', 'w/MS')
title(['Mag of echo # ' num2str(necho) ', shot # ' num2str(nshot)], 'fontsize', 18)


