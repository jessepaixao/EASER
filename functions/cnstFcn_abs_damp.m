function [c,ceq] = cnstFcn_abs_damp(x,ms,ks,ma,omega)

    ka=x(1);
    ca=x(2);
   
    Nm=length(ks);
    H=zeros(Nm,length(omega));
    
    for i=1:Nm
        ks_rnd=ks(i);
        for j=1:length(omega)

            lamb=omega(j)/(sqrt(ks_rnd/ms));
            gamma=sqrt(ka/ma)/sqrt(ks_rnd/ms);
            mu=ma/ms;
            eta=ca/(2*sqrt(ma*ka));

            % Receptance of coupled system
            num=(gamma^2-lamb^2)^2+(2*eta*gamma*lamb)^2;
            den=(((1-lamb^2)*(gamma^2-lamb^2)-mu*gamma^2*lamb^2)^2+(2*eta*gamma*lamb)^2*(1-mu*lamb^2-lamb^2)^2);
            H(i,j)=sqrt(num/den);


        end
    end
    
    h_max=max(H,[],1);
    
    f=norm(h_max,Inf);

    c = f-x(3);
    ceq = [ ];

end