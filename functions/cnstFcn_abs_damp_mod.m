function [c,ceq] = cnstFcn_abs_damp_mod(x,m1,k1,c1,m2,omega,c2)

    k2=x(1);
%     c2=x(2);
   
    Nm=length(k1);
    Hc=zeros(Nm,length(omega));
    h_max=[];
    for i=1:Nm
        ks_rnd=k1(i);
        for k=1:length(omega)
             dv2=c2;
             dv1=c1;
    
             Hc(i,k)=abs((-omega(k)^2*m2+k2+1j*dv2)/((-omega(k)^2*m1+1j*dv1+1j*dv2+ks_rnd+k2)*(-omega(k)^2*m2+k2+1j*dv2)-(1j*dv2+k2).^2));  
            
%             lamb=omega(j)/(sqrt(ks_rnd/ms));
%             gamma=sqrt(ka/ma)/sqrt(ks_rnd/ms);
%             mu=ma/ms;
%             eta=ca/(2*sqrt(ma*ka));
% 
%             % Receptance of coupled system
%             num=(gamma^2-lamb^2)^2+(2*eta*gamma*lamb)^2;
%             den=(((1-lamb^2)*(gamma^2-lamb^2)-mu*gamma^2*lamb^2)^2+(2*eta*gamma*lamb)^2*(1-mu*lamb^2-lamb^2)^2);
%             H(i,j)=sqrt(num/den);


        end
        
%         norm(abs(Hc(i,:)),Inf)
%         h_max(i)=norm(Hc(i,:),Inf)*1E4;
    end
    
    h_max=max(Hc,[],1)*1e4;
%     h_max
    f=max(h_max);

    c = f-x(2);
    ceq = [ ];

end