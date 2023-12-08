%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of using Abers and Hacker 2016 with the VBRc:
% - unrelaxed moduli and density calculated with Abers and Hacker for a rock
%   composition of 90% Forsterite, 10% Fayalite (by volume)
% - relaxed properties calculated with the VBRc
% - based on `rockvelcalculate.m` from Abers and Hacker
%
% Setup: 
% - installation of the vbrc, path set with `vbrdir` environment variable 
% - code from Abers and Hacker2016: download and set the `ABERSHACKER16` environment
%   variable to the path of the unzipped directory. Abers and Hacker code available in
%   the supplement of:
%
%   Abers, G. A., and Hacker, B. R. (2016), A MATLAB toolbox and Excel workbook for
%   calculating the densities, seismic wave speeds, and major element composition of
%   minerals and rocks at pressure and temperature, Geochem. Geophys. Geosyst., 17,
%   616–624, doi:10.1002/2015GC006171.
%   https://doi.org/10.1002/2015GC006171
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
addpath(genpath(getenv("ABERSHACKER16")))
addpath(getenv('vbrdir'))
vbr_init();

% Set notional rock composition
mins={'fo', 'fa'};
modes=[90 10];   % These are volume fraction modes

% Set conditions for anharmonic calculation
nT = 10;
nP = 15;
T_K1d = linspace(1000, 1773, nT);
P_GPa1d = linspace(1, 3, nP);

[T_K, P_GPa] = meshgrid(T_K1d, P_GPa1d);

% Load mineral database
[minpropar, compar]=ah16_loaddb('AbersHackerMacroJan2016.txt');
if isempty(minpropar)
    disp('ERROR on AH16_LOADDB')
    return
end

% initialize output arrays for moduli and density
G = zeros(size(T_K));
K = zeros(size(T_K));
rho = zeros(size(T_K));

for i_T = 1:nT
    for i_P = 1:nP
        tt = T_K(i_P, i_T)-273.;  %deg C
        pp = P_GPa(i_P, i_T);  % GPa

        %% Calcuate - single P,T  -- also get rhos
        [modu,modhsm,modhsp,modvrm,modvrp,rhos]=ah16_rockvel(tt,pp, minpropar, mins,modes);

        % first calculate correct molar Fo90 volume fractions: want 90 mol% fo
        %    for high accuracy this requires 2 calls to ah16_rockvel, the first
        %    to get rhos for nominal modes (above), then a second with correct vol% modes.
        k_fo=find(strcmp('fo',{minpropar.name}));
        k_fa=find(strcmp('fa',{minpropar.name}));
        modecor=[minpropar(k_fo).gfw minpropar(k_fa).gfw]./[rhos(k_fo) rhos(k_fa)]; % cm3/mol per phase

        modes_c = modes.*modecor;
        [modu,modhsm,modhsp,modvrm,modvrp,compos]=ah16_rockvel(tt,pp, minpropar, mins,modes_c,compar);

        G(i_P, i_T) = modu.g*1e9;
        K(i_P, i_T) = modu.k*1e9;
        rho(i_P, i_T) = modu.r;
    end
end

% setup the VBRc calculation for anelastic properties

VBR.in.elastic.methods_list={'anharmonic';};
VBR.in.anelastic.methods_list={'eburgers_psp';'andrade_psp';};

% set state variables: first, copy over Abers and Hacker inputs and outputs
VBR.in.SV.rho = rho;
VBR.in.SV.P_GPa = P_GPa; % pressure [GPa]
VBR.in.SV.T_K = T_K; % temperature [K]
VBR.in.elastic.Gu_TP = G;
VBR.in.elastic.Ku_TP = K;

% and the remaining state variables
VBR.in.SV.f = logspace(-5,-1, 50); % frequency [Hz]
VBR.in.SV.sig_MPa = .1 * ones(size(T_K)); % differential stress [MPa]
VBR.in.SV.phi = zeros(size(T_K)); % melt fraction
VBR.in.SV.dg_um = 0.01 * 1e6 * ones(size(T_K)); % grain size [um]

% and the calculation
VBR = VBR_spine(VBR);

% plot log10(Q) for a single frequency
figure('PaperPosition',[0,0,14,6],'PaperPositionMode','manual')
subplot(1,2,1)
contourf(T_K1d, P_GPa1d, log10(squeeze(VBR.out.anelastic.eburgers_psp.Q(:,:,end))), 20)
caxis([1 7])
colorbar()
xlabel("Temperature [K]")
ylabel("Pressure [GPa]")
title("log10(Q), extended burgers pseudo-period (JF10)")
cmapname = "cubehelix";
colormap(cmapname)
subplot(1,2,2)
contourf(T_K1d, P_GPa1d, log10(squeeze(VBR.out.anelastic.andrade_psp.Q(:,:,end))), 20)
caxis([1 7])
colorbar()
colormap(cmapname)
xlabel("Temperature [K]")
ylabel("Pressure [GPa]")
title("log10(Q), andrade pseudo-period (JF10)")
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(gcf,'./AbersHackerVBRc_Q.png')

% plot frequency dependence at 1500 C, 3 GPa
figure('PaperPosition',[0,0,4,4],'PaperPositionMode','manual')
iT = nT;
iP = nP;
f = VBR.in.SV.f;
semilogx(f, squeeze(VBR.out.anelastic.andrade_psp.M(iP, iT, :)/1e9), 'k','displayname', 'andrade', 'linewidth',1.5)
hold all
semilogx(f, squeeze(VBR.out.anelastic.eburgers_psp.M(iP, iT, :)/1e9),'--r', 'displayname', 'eburgers', 'linewidth',1.5)
legend('location', 'NorthWest')
xlabel('frequency [Hz]')
ylabel('M [GPa]')
title('1500 C, 3 GPa')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(gcf,'./AbersHackerVBRc_M_v_f.png')


