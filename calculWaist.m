% by yann leventoux 2025
close all
clear all

fiberGUI

function fiberGUI
    % Main interface
    fig = uifigure('Name','Fundamental mode calculator','Position',[100 100 420 380]);

    % Input fields
    uilabel(fig,'Text','Core diameter (µm):','Position',[20 310 160 22]);
    dCoreField = uieditfield(fig,'numeric','Position',[200 310 100 22],'Value',8);

    uilabel(fig,'Text','Wavelength (µm):','Position',[20 270 160 22]);
    lambdaField = uieditfield(fig,'numeric','Position',[200 270 100 22],'Value',1.55);

    uilabel(fig,'Text','Given parameter:','Position',[20 230 160 22]);
    typeMenu = uidropdown(fig,'Items',{'NA','delta n','n1'}, ...
                          'Position',[200 230 100 22]);

    uilabel(fig,'Text','Value:','Position',[20 190 160 22]);
    valueField = uieditfield(fig,'numeric','Position',[200 190 100 22]);

    % Default values
    typeMenu.ValueChangedFcn = @(src,event) setDefaultValue(src,valueField);
    typeMenu.Value = 'NA';
    valueField.Value = 0.1;

    % Calculate button
    uibutton(fig,'Text','Calculate','Position',[160 150 100 30],...
        'ButtonPushedFcn',@(btn,event) calculate(dCoreField,lambdaField,typeMenu,valueField,fig));

    % Results box
    uitextarea(fig,'Editable','off','Position',[20 20 380 120],'Value',{'Results:'});
end

function setDefaultValue(typeMenu,valueField)
    switch typeMenu.Value
        case 'NA'
            valueField.Value = 0.1;
        case 'delta n'
            valueField.Value = 1e-3;
        case 'n1'
            valueField.Value = 1.45;
    end
end

function calculate(dCoreField,lambdaField,typeMenu,valueField,fig)
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
    rcl = 10*a;
    [MFD_gauss, MFD_pet, Aeff_rig, neff, rp, Ep] = computeMode(a,rcl,n1,n2,lambda);

    % Approx effective area from Gaussian
    Aeff_approx = pi*(MFD_gauss/2)^2;

    % V-number
    V = 2*pi*a*NA/lambda;

    % Display results
    txt = sprintf(['n1 = %.6f, n2 = %.6f\n' ...
                   'neff = %.6f\n' ...
                   'MFD (Gaussian 1/e²) = %.3f µm\n' ...
                   'MFD (Petermann II)  = %.3f µm\n' ...
                   'Aeff (pi*w0²)  = %.3f µm²\n' ...
                   'Aeff (field integral)  = %.3f µm²\n' ...
                   'V-number = %.3f'], ...
                   n1,n2,neff,MFD_gauss*1e6,MFD_pet*1e6,Aeff_approx*1e12,Aeff_rig*1e12,V);
    res = findobj(fig,'Type','uitextarea');
    res.Value = {'Results:',txt};

    % ---- Plot field and MFDs ----
    figure(1);
    clf
    plot(rp*1e6, abs(Ep).^2/max(abs(Ep).^2),'b','LineWidth',1.5); hold on;
    xline(MFD_gauss/2*1e6,'--r','W0 (Gaussian)');
    xline(MFD_pet/2*1e6,'--g','W0 (Petermann)');
    xlabel('Radius (µm)');
    ylabel('Normalized intensity');
    title('LP01 mode profile and MFDs');
    legend('Mode intensity','Gaussian MFD radius','Petermann MFD radius');
    grid on;
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

function [MFD_gauss, MFD_pet, Aeff, neff, rp, Ep] = computeMode(a,rcl,n1,n2,lambda)
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
        U1 = besselLP01(r,n1,n2,rco,nGuess1,lambda);
        U1p = gradient(gradient(U1,r),r);
        U2 = besselLP01(r,n1,n2,rco,nGuess2,lambda); 
        U2p = gradient(gradient(U2,r),r);
        if min(U1p)>min(U2p)
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

    % --- Petermann MFD
    num = trapz(rp, (rp.^2) .* abs(Ep).^2 );
    num = num^2;
    den = trapz(rp, (rp.^3) .* abs(Ep).^4 );
    wp = sqrt(2 * num / den);
    MFD_pet = 2*wp;

    % --- Rigorous effective area
    num2 = trapz(rp, abs(Ep).^2 .* rp);
    num2 = num2^2;
    den2 = trapz(rp, abs(Ep).^4 .* rp);
    Aeff = 2*pi * num2/den2;
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
