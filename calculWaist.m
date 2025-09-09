% by Yann Leventoux 2025
close all
clear all

fibreGUI

function fibreGUI
    % Interface principale
    fig = uifigure('Name','Calcul du mode LP01','Position',[100 100 420 340]);

    % Champs de saisie
    uilabel(fig,'Text','Diamètre du coeur (µm):','Position',[20 270 160 22]);
    dCoreField = uieditfield(fig,'numeric','Position',[200 270 100 22],'Value',8);

    uilabel(fig,'Text','Longueur d''onde (µm):','Position',[20 230 160 22]);
    lambdaField = uieditfield(fig,'numeric','Position',[200 230 100 22],'Value',1.55);

    uilabel(fig,'Text','Paramètre donné :','Position',[20 190 160 22]);
    typeMenu = uidropdown(fig,'Items',{'NA',[char(916) 'n'],'n1'}, ...
                          'Position',[200 190 100 22]);

    uilabel(fig,'Text','Valeur :','Position',[20 150 160 22]);
    valueField = uieditfield(fig,'numeric','Position',[200 150 100 22]);

    % Callback pour mettre valeur par défaut
    typeMenu.ValueChangedFcn = @(src,event) setDefaultValue(src,valueField);

    % Initialisation : NA = 0.1
    typeMenu.Value = 'NA';
    valueField.Value = 0.1;

    % Bouton de calcul
    uibutton(fig,'Text','Calculer','Position',[160 90 100 30],...
        'ButtonPushedFcn',@(btn,event) calculer(dCoreField,lambdaField,typeMenu,valueField,fig));

    % Zone de résultat
    uitextarea(fig,'Editable','off','Position',[20 20 380 60],'Value',{'Résultats :'});
end

function setDefaultValue(typeMenu,valueField)
    switch typeMenu.Value
        case 'NA'
            valueField.Value = 0.1;
        case [char(916) 'n']   % ?n
            valueField.Value = 1e-3;
        case 'n1'
            valueField.Value = 1.45;
    end
end

function calculer(dCoreField,lambdaField,typeMenu,valueField,fig)
    % Récupération paramètres
    a = dCoreField.Value/2*1e-6;   % rayon coeur [m]
    lambda = lambdaField.Value*1e-6; % [m]

    % Indice de la gaine via Sellmeier (silice pure)
    n2 = indiceSiGe(lambda,0); 

    % Détermination de n1
    val = valueField.Value;
    switch typeMenu.Value
        case 'NA'
            NA = val;
            n1 = sqrt(n2^2 + NA^2);
        case [char(916) 'n']
            dn = val;
            n1 = n2 + dn;
        case 'n1'
            n1 = val;
    end

    % NA utile (toujours bien d'avoir)
    NA = sqrt(n1^2 - n2^2);

    % Calcul waist
    rcl = 3*a;
    W0 = calculW0(a,rcl,n1,n2,lambda);

    % Aire effective approx gaussienne
    Aeff = pi*W0.^2;

    % V-number
    V = 2*pi*a*NA/lambda;

    % Résultats
    txt = sprintf(['n1 = %.6f, n2 = %.6f\n' ...
                   'W0 = %.3f µm   Aeff = %.3f µm²\n' ...
                   'V = %.3f'], ...
                   n1,n2,W0*1e6,Aeff*1e12,V);
    res = findobj(fig,'Type','uitextarea');
    res.Value = {'Résultats :',txt};
end


function indice = indiceSiGe( lambda, XGe )
    GeB1 = 0.80686642;
    GeB2 = 0.71815848;
    GeB3 = 0.85416831;
    GeC1 = ((0.68972606e-1)^2)*1e-12;
    GeC2 = (0.15396605^2)*1e-12;
    GeC3 = ((0.11841931e2)^2)*1e-12;
    SiB1 = 0.696166300;
    SiB2 = 0.407942600;
    SiB3 = 0.897479400;
    SiC1 = 4.67914826e-3*1e-12;
    SiC2 = 1.35120631e-2*1e-12;
    SiC3 = 97.9340025*1e-12;
    %les A du papier c'est mes B et les Q² du papier c'est mes C
    indice = sqrt(1 + ((SiB1+XGe*(GeB1-SiB1))*lambda^2)/(lambda^2-(sqrt(SiC1)+XGe*(sqrt(GeC1)-sqrt(SiC1)))^2)...
        +((SiB2+XGe*(GeB2-SiB2))*lambda^2)/(lambda^2-(sqrt(SiC2)+XGe*(sqrt(GeC2)-sqrt(SiC2)))^2)...
        + ((SiB3+XGe*(GeB3-SiB3))*lambda^2)/(lambda^2-(sqrt(SiC3)+XGe*(sqrt(GeC3)-sqrt(SiC3)))^2));
end
function W0 = calculW0(a,rcl,n1,n2,lambda)
    nbTry = 30; %nb try pour dichotomie
    wl = lambda;
    rco = a;
    r = linspace(-rcl,rcl,30000);
    neff1 = n1;
    neff2 = n2;
    neff3 = (2*neff1+neff2)./3;
    neff4 = (neff1+2*neff2)./3;
    for nb = 1:nbTry
        U1 = besselLP01(r,n1,n2,rco,neff1,wl);
        U2 = besselLP01(r,n1,n2,rco,neff2,wl);
        U3 = besselLP01(r,n1,n2,rco,neff3,wl);
        U4 = besselLP01(r,n1,n2,rco,neff4,wl);  
        temp = [max(diff(U1)) max(diff(U2)) max(diff(U3)) max(diff(U4))];
        nmin = find(temp == min(temp));
        if nmin == 1
            neff1 = neff1;
            neff2 = neff3;
            neff3 = (2*neff1+neff2)./3;
            neff4 = (neff1+2*neff2)./3;
        elseif nmin == 2
            neff1 = neff4;
            neff2 = neff2;
            neff3 = (2*neff1+neff2)./3;
            neff4 = (neff1+2*neff2)./3; 
        elseif nmin == 3
            neff1 = neff1;
            neff2 = neff4;
            neff3 = (2*neff1+neff2)./3;
            neff4 = (neff1+2*neff2)./3;        
        else
            neff1 = neff3;
            neff2 = neff2;
            neff3 = (2*neff1+neff2)./3;
            neff4 = (neff1+2*neff2)./3;
        end
    end
    U = besselLP01(r,n1,n2,rco,neff1,wl);
    I = abs(U).^2;
    I = I./max(I);
    temp = r(find(I>0.1353));
    neff = neff1;
    W0 = max(temp);
end

function U = besselLP01(r,n1,n2,rco,neff,wl)
    ko = 2*pi/wl;
    beta = 2*pi*neff/wl;
    kT = sqrt(n1^2*ko^2-beta^2);
    gamma = sqrt(beta^2-n2^2*ko^2);
    coefMerde = besselj(0,kT.*rco)/besselk(0,gamma*rco);
    U = r.*0;
    U(find(abs(r)<rco)) = besselj(0,kT.*r(find(abs(r)<rco)));
    U(find(abs(r)>rco)) = coefMerde.*real(besselk(0,gamma.*r(find(abs(r)>rco))));
end
