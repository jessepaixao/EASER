function h=recept2MCK_opt(omega,m1,k1,c1,m2,k2,c2)

    
    
    for i=1:length(omega)

        dv2=c2;
        dv1=c1;
    
        Hc(i)=((-omega(i)^2*m2+k2+1j*dv2)/((-omega(i)^2*m1+1j*dv1+1j*dv2+k1+k2)*(-omega(i)^2*m2+k2+1j*dv2)-(1j*dv2+k2).^2));  
    

    end
    
    h=norm(Hc,Inf);
    
    
