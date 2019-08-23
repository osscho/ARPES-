function ret = dataslice(kx,ky,Ekin,V, angles_x, angles_y, k_surf_x, k_surf_y, E_surf)

    %addpath('./data');


    
    x_steps = size(k_surf_x,2);
    y_steps = size(k_surf_y,1);

    xi=zeros(y_steps,x_steps);
    yi=xi;
    zi=xi;

    maxE=max(max(max(Ekin)));
    minE=min(min(min(Ekin)));
    E_step=(maxE-minE)/(size(kx,3)-1);

    maxax=angles_x(size(angles_x,2));
    minax=angles_x(1);
    ax_step=(maxax-minax)/(size(angles_x,2)-1);

    maxay=angles_y(size(angles_y,2));
    minay=angles_y(1);
    ay_step=(maxay-minay)/(size(angles_y,2)-1);
    %close all;
    %hold on;
    for itx = 1:y_steps
        for ity = 1:x_steps
            a_x =360/2/pi* asin(k_surf_x(itx,ity)/(0.5123*sqrt(E_surf(itx,ity))));
            a_y =360/2/pi* asin(k_surf_y(itx,ity)/(0.5123*sqrt(E_surf(itx,ity))));
            %plot(k_surf_y(itx,ity))
            zi(itx,ity)=(E_surf(itx,ity)-minE)/E_step;
            yi(itx,ity)=(a_x-minax)/ax_step;
            xi(itx,ity)=(a_y-minay)/ay_step;
            if not(xi(itx,ity) == conj(xi(itx,ity)))
                xi(itx,ity)=nan;
            end
            if not(yi(itx,ity) == conj(yi(itx,ity)))
                yi(itx,ity)=nan;
            end
            if not(zi(itx,ity) == conj(zi(itx,ity)))
                zi(itx,ity)=nan;
            end
               
        end
    end
    %hold off;
    
    ret=interp3(V,xi,yi,zi,'cubic');
end

