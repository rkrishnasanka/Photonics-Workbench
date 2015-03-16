function plotAWGTEMode(beta_field,k,ns,nf,nc,h_s,h_f,h_c)

    %plot field profiles

    A_s = 1;

    gamma_s = sqrt(beta_field^2-(ns*k)^2);
    gamma_c = sqrt(beta_field^2-(nc*k)^2);
    kappa_f = sqrt((nf*k)^2-beta_field^2);

    NN = 100;

    thickness = h_f + h_c + h_s;

    step = thickness/NN;

    x_step = 0;

    for i = 1:NN+1

        x_step(i+1) = x_step(i) + step;

        if (x_step(i) <= h_s);

            E_y(i) = A_s*exp(gamma_s*(x_step(i) - h_s));

        elseif (h_s <= x_step(i)) && (x_step(i) <= h_s + h_f);

            E_y(i) = A_s*(cos(kappa_f*(x_step(i) - h_s)) + gamma_s*sin(kappa_f*(x_step(i)- h_s))./kappa_f);

        elseif (h_s + h_f <= x_step(i)) && (x_step(i) <= thickness);

            %The more complicated E_y(i) is formed here because of the boundry
            %condition matching we have to do to make Amplitudes for the fields
            %have to be matched
            E_y(i) = A_s*(cos(kappa_f*h_f) + gamma_s * sin(kappa_f*h_f)/kappa_f)*exp(-gamma_c*(x_step(i) - h_f - h_s)); 

        end
    end

    % 

    x = 0.0:step:thickness; %X Axis coordinates of plot points

    h=plot(x,E_y);
    vline(h_s,'r','s-f');
    vline(h_s + h_f,'r','f-c');
    xlabel('x (microns)','FontSize',22);
    ylabel('TE electric field','FontSize',22);
    set(h,'LineWidth',1.5); %new thickness of plotting lines
    set(gca,'FontSize',22); %new size of tick marks on both axes
    grid on

end