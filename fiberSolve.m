%calcul of neff, disp, vg, MFD, Aeff of fondamental mode via bessel
%resolution
%by yann Leventoux 2025
close all
clear all

fiberGUI

function fiberGUI
    % Main interface
    fig = uifigure('Name','Fundamental mode calculator','Position',[100 120 420 420]);

    % Input fields
    uilabel(fig,'Text','Core diameter (µm):','Position',[20 360 160 22]);
    dCoreField = uieditfield(fig,'numeric','Position',[200 360 100 22],'Value',8.2);

    uilabel(fig,'Text','Wavelength (µm):','Position',[20 320 160 22]);
    lambdaField = uieditfield(fig,'numeric','Position',[200 320 100 22],'Value',1.55);

    uilabel(fig,'Text','Given parameter:','Position',[20 280 160 22]);
    typeMenu = uidropdown(fig,'Items',{'NA','delta n','n1'}, ...
                          'Position',[200 280 100 22]);

    uilabel(fig,'Text','Value:','Position',[20 2140 160 22]);
    valueField = uieditfield(fig,'numeric','Position',[200 240 100 22]);

    % Default values
    typeMenu.ValueChangedFcn = @(src,event) setDefaultValue(src,valueField);
    typeMenu.Value = 'NA';
    valueField.Value = 0.14;

    % Checkbox pour calculer la dispersion
    dispersionBox = uicheckbox(fig, ...
        'Text','neff(lambda), Vg, dispersion', ...
        'Position',[20 200 200 22], ...
        'Value',false);

    % Calculate button
    uibutton(fig,'Text','Calculate','Position',[160 160 100 30],...
        'ButtonPushedFcn',@(btn,event) calculate(dCoreField,lambdaField,typeMenu,valueField,dispersionBox,fig));

    % Results box (agrandie)
    uitextarea(fig,'Editable','off','Position',[20 10 380 140],'Value',{'Results:'});
end

function setDefaultValue(typeMenu,valueField)
    switch typeMenu.Value
        case 'NA'
            valueField.Value = 0.14;
        case 'delta n'
            valueField.Value = 1e-3;
        case 'n1'
            valueField.Value = 1.45;
    end
end

function calculate(dCoreField,lambdaField,typeMenu,valueField,dispersionBox,fig)
    % Parameters
    a = dCoreField.Value/2*1e-6;      % core radius [m]
    lambda = lambdaField.Value*1e-6;  % wavelength [m]

    % Cladding index (pure silica via Sellmeier)
    n2 = silicaIndex(lambda,0); 

    % Core index determination
    val = valueField.Value;
    switch typeMenu.Value
        case 'NA'
            NA = val;
            n1 = sqrt(n2^2 + NA^2);
        case 'delta n'
            dn = val;
            n1 = n2 + dn;
        case 'n1'
            n1 = val;
    end

    % Numerical aperture
    NA = sqrt(n1^2 - n2^2);

    % Mode calculation
    rcl = 5*a;
    [MFD_gauss, MFD_pet, MFD_4sigma, Aeff_rig, neff, rp, Ep] = computeMode(a,rcl,n1,n2,lambda);

    % Approx effective area from Gaussian
    Aeff_approx = pi*(MFD_gauss/2)^2;

    % V-number
    V = 2*pi*a*NA/lambda;

    % ---- Résultats texte ----
    txt = sprintf(['n1 = %.6f, n2 = %.6f\n' ...
                   'neff = %.6f\n' ...
                   'MFD (Gaussian 1/e²) = %.3f µm\n' ...
                   'MFD (Petermann 2)  = %.3f µm\n' ...
                   'MFD (near-field rms) = %.3f µm\n' ...
                   'Aeff (pi*w0²)  = %.3f µm²\n' ...
                   'Aeff (field integral)  = %.3f µm²\n' ...
                   'V-number = %.3f'], ...
                   n1,n2,neff,MFD_gauss*1e6,MFD_pet*1e6,MFD_4sigma*1e6,...
                   Aeff_approx*1e12,Aeff_rig*1e12,V);

    % ---- Optionnel : calcul dispersion ----
    if dispersionBox.Value
        D_interp = calcDispersion(typeMenu.Value,valueField.Value,a,lambda);
        txt = sprintf('%s\nDispersion at lambda=%.2f µm : %.2f ps/nm/km', ...
                      txt, lambda*1e6, D_interp);
    end

    res = findobj(fig,'Type','uitextarea');
    res.Value = {'Results:',txt};

    % ---- Plot field and MFDs ----
    figure(1);clf
    plot(rp*1e6, abs(Ep).^2/max(abs(Ep).^2),'b','LineWidth',1.5); hold on;
    xline(MFD_gauss/2*1e6,'--r','W0 (Gaussian)');
    xline(MFD_pet/2*1e6,'--g','W0 (Petermann)');
    xline(MFD_4sigma/2*1e6,'--m','W0 (4 sigma))');
    xlabel('Radius (µm)','FontSize',12);
    ylabel('Normalized intensity','FontSize',12);
    title('LP01 mode profile and MFDs','FontSize',13);
    legend('Mode intensity','Gaussian MFD radius','Petermann MFD radius','4\sigma MFD radius');
    grid on; box on;
end

function index = silicaIndex(lambda, XGe)
    % Sellmeier equation for pure silica / Ge-doped silica
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
    index = sqrt(1 + ((SiB1+XGe*(GeB1-SiB1))*lambda^2)/(lambda^2-(sqrt(SiC1)+XGe*(sqrt(GeC1)-sqrt(SiC1)))^2)...
        +((SiB2+XGe*(GeB2-SiB2))*lambda^2)/(lambda^2-(sqrt(SiC2)+XGe*(sqrt(GeC2)-sqrt(SiC2)))^2)...
        + ((SiB3+XGe*(GeB3-SiB3))*lambda^2)/(lambda^2-(sqrt(SiC3)+XGe*(sqrt(GeC3)-sqrt(SiC3)))^2));
end

function [MFD_gauss, MFD_pet, MFD_4sigma, Aeff, neff, rp, Ep] = computeMode(a,rcl,n1,n2,lambda)
    rco = a;
    r = linspace(-rcl,rcl,3000);
    Nguess = 100;
    
    % Estimate neff_min from cutoff condition V=2.405
    neffVect = linspace(n2,n1,10000);
    Veff = 2*pi/lambda*rco.*sqrt(n1^2-neffVect.^2);
    [~, idx] = min(abs(Veff-2.405));
    neffMin = neffVect(idx);
   
    a = n1;
    b = neffMin;
    for nb = 1:Nguess
        nGuess1 = a+(b-a)/3;
        nGuess2 = b-(b-a)/3;
        eq1 = eqDispCalulation(n1,n2,rco,nGuess1,lambda);
        eq2 = eqDispCalulation(n1,n2,rco,nGuess2,lambda);
        if eq1*eq2<0
            b = nGuess2;
        else
            a = nGuess1;
        end
    end

    neff = nGuess1;
    U = besselLP01(r,n1,n2,rco,neff,lambda);
    rp = r(r>=0);
    Ep = U(r>=0);

    % Normalize energy
    normFactor = sqrt(trapz(rp, abs(Ep).^2 .* rp) * 2*pi);
    Ep = Ep ./ normFactor;

    % --- Gaussian MFD (1/e²)
    I = abs(Ep).^2;
    I = I./max(I);
    temp = rp(I>0.1353);
    W0 = max(temp);
    MFD_gauss = 2*W0;

    % --- Petermann 2 MFD
    dEp = gradient(Ep,rp);
    num = trapz(rp, rp .* Ep.^2 );
    den = trapz(rp, rp .* dEp.^2 );
    wp = sqrt(2 * num / den);
    MFD_pet = 2*wp;

    % --- 4sigma MFD or near-field rms
    num = trapz(rp, (rp.^3) .* abs(Ep).^2);
    den = trapz(rp, (rp) .* abs(Ep).^2);
    wsig = 2*sqrt(num/den);
    MFD_4sigma = 2*wsig/sqrt(2);

    % --- Rigorous effective area
    num = 2*pi * (trapz(rp, (Ep.^2) .* rp))^2;
    den = trapz(rp, (Ep.^4) .* rp);
    Aeff = num / den;
end

function U = besselLP01(r,n1,n2,rco,neff,wl)
    ko = 2*pi/wl;
    beta = 2*pi*neff/wl;
    kT = sqrt(n1^2*ko^2-beta^2);
    gamma = sqrt(beta^2-n2^2*ko^2);
    coef = besselj(0,kT.*rco)/besselk(0,gamma*rco);
    U = zeros(size(r));
    U(abs(r)<rco) = besselj(0,kT.*r(abs(r)<rco));
    U(abs(r)>=rco) = coef.*real(besselk(0,gamma.*abs(r(abs(r)>=rco))));
end

function D_interp = calcDispersion(paramType,paramValue,a,lambda0)
    % Paramètres spectre
    lambdaMin = 0.55e-6;
    lambdaMax = 2.4e-6;
    Nlambda = 100;
    lambdaVec = linspace(lambdaMin,lambdaMax,Nlambda);

    % Constantes
    c = 3e8;
    Nguess = 100;
    rco = a;
    neff = zeros(1,Nlambda);

    % --- Calcul de neff(?)
    for nb = 1:Nlambda
        lambda_temp = lambdaVec(nb);
        n2 = silicaIndex(lambda_temp,0); % cladding
        switch paramType
            case 'delta n'
                n1 = n2 + paramValue;
            case 'NA'
                NA = paramValue;
                n1 = sqrt(n2^2 + NA^2);
            case 'n1'
                n1l0 = paramValue;
                n2l0 = silicaIndex(lambda0,0);
                dn =n1l0 - n2l0;
                n1 = n2 + dn;       
        end
        neff(nb) = dichotomieBessel(n1,n2,rco,lambda_temp,Nguess);
    end

    % --- Dérivées numériques
    dneff_dlambda = gradient(neff, lambdaVec);     
    d2neff_dlambda2 = gradient(dneff_dlambda, lambdaVec);

    % --- Vitesse de groupe
    ng = neff - lambdaVec .* dneff_dlambda; % indice de groupe
    vG = c ./ ng; 

    % --- Dispersion (ps/nm/km)
    D = -(lambdaVec ./ c) .* d2neff_dlambda2; 
    D_psnmkm = D * 1e6;

    % --- Plot neff(?) ---
    figure(2);clf
    plot(lambdaVec*1e6, neff,'-b','LineWidth',1.5);
    xlabel('Wavelength \lambda (µm)','FontSize',12);
    ylabel('n_{eff}','FontSize',12);
    title('Effective index vs wavelength','FontSize',13);
    grid on; box on;

    % --- Plot vG(?) (sans 1er et dernier point) ---
    figure(3);clf
    plot(lambdaVec(2:end-1)*1e6, vG(2:end-1)/1e8,'-g','LineWidth',1.5);
    xlabel('Wavelength \lambda (µm)','FontSize',12);
    ylabel('v_g (10^8 m/s)','FontSize',12);
    title('Group velocity vs wavelength','FontSize',13);
    grid on; box on;

    % --- Plot D(?) (sans 2 premiers et 2 derniers points) ---
    figure(4);clf
    plot(lambdaVec(3:end-2)*1e6, D_psnmkm(3:end-2),'-r','LineWidth',1.5);
    xlabel('Wavelength \lambda (µm)','FontSize',12);
    ylabel('D (ps/nm/km)','FontSize',12);
    title('Chromatic dispersion vs wavelength','FontSize',13);
    grid on; box on;

    % --- Interpolation à lambda0 ---
    D_interp = interp1(lambdaVec, D_psnmkm, lambda0,'linear','extrap');
end

function neff = dichotomieBessel(n1,n2,rco,lambda,Nguess)
    a = n1;  
    % Estimate neff_min from cutoff condition V=2.405
    neffVect = linspace(n2,n1,10000);
    Veff = 2*pi/lambda*rco.*sqrt(n1^2-neffVect.^2);
    [~, idx] = min(abs(Veff-2.405));
    neffMin = neffVect(idx);  
    b = neffMin;
    for nb = 1:Nguess
        nGuess1 = a+(b-a)/3;
        nGuess2 = b-(b-a)/3;
        eq1 = eqDispCalulation(n1,n2,rco,nGuess1,lambda);
        eq2 = eqDispCalulation(n1,n2,rco,nGuess2,lambda);
        if eq1*eq2<0
            b = nGuess2;
        else
            a = nGuess1;
        end
    end
    neff = nGuess1;
end

function eq = eqDispCalulation(n1,n2,rco,neff,lambda)
    k0 = 2*pi/lambda;
    beta = k0*neff;
    U = rco .* sqrt((k0*n1).^2 - beta.^2);
    W = rco .* sqrt(beta.^2 - (k0*n2).^2);  
    lhs = besselj(0,U) ./ (U .* besselj(1,U));   % côté cœur
    rhs = (besselk(0,W)) ./ ((W .* besselk(1,W)));   % côté gaine 
    eq = lhs - rhs;   % doit s’annuler aux racines
end

