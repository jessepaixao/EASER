function h=recept2MCK(omega,m1,k1,c1,m2,k2,c2)

    N=length(omega);

    % Initialization of receptances
    h=zeros(1,N);

    for i=1:length(omega)

        
        dv2=c2;
        dv1=c1;
    
        h(i)=((-omega(i)^2*m2+k2+1j*dv2)/((-omega(i)^2*m1+1j*dv1+1j*dv2+k1+k2)*(-omega(i)^2*m2+k2+1j*dv2)-(1j*dv2+k2).^2));  
    

    end
    
    
