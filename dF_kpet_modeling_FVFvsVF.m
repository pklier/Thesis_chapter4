%models dF/F of Voltage Fluors as a function of membrane potential based on measured and
%calculated PeT parameters

hold on
%clear all;

dEox = 0.129;
dEred = -2.02;
dEoo = 2.38;
dHOMO = -0.231;
%setting up constants for energy of electron transfer equation; can be modified to use
%HOMO-HOMO or LUMO-LUMO gap. Currently set to Evan's values for VF2.1Cl
%(note that redox potentials are measured in volts but since this is a
%transit of one electron this is directly equivalent to electron volts in
%energy)

r = 2.2;
%aniline-chromophore distance (nm)
Vmem = [-300:0.1:300];
%sets up a vector to calculate kpet for a wide range of membrane potentials

%theta = 35.3;
theta = [0:2.5:90];
weight = [0.0036    0.0101    0.0165    0.0232    0.0299    0.0362    0.0417    0.0442    0.0480    0.0524    0.0544 0.0553    0.0554    0.0535    0.0522    0.0499    0.0491    0.0466    0.0423    0.0371    0.0320    0.0282 0.0239    0.0203    0.0180    0.0162    0.0138    0.0118    0.0097    0.0074    0.0058    0.0042    0.0033 0.0021    0.0011    0.0006         0];

%angle between wire and electric field; should be modeled computationally
%or this program can be modified to calculate this from measured dF/F.
%Currently set to the computationally derived distribution of angles for
%VF2.1Cl from Rishi's paper
dmem = 4;
%thickness of the membrane in question (nm); 4 nm is a ballpark number. Not
%dubelco's modified eagle media
w = (r.*Vmem./(1000*dmem)).*sum(weight.*cosd(theta));
%w = (r.*Vmem./(1000*dmem)).*cosd(theta);

%calculates work to move electron in PeT (signs and charge accounted for, should be net positive so watch out
%for cosine range) based on above parameters (takes cosine in degrees);
%1000 in denominator accounts for mV to V conversion
%sanity check: work is positive when Vmem is positive (takes energy to move
%an electron from positive to negative) and work is negative when Vmem is
%negative (electron wants to go from negative to positive)

dGpet = dEox - dEred - dEoo + w;
%dGpet = dHOMO +w;
%Rehm-Weller equation to calculate delta G of PeT as a function of Vmem; we will now use this to
%calculate kpet as a function of Vmem which we will then use to calculate
%dF/F as a function of Vmem
%+w is because of current sign convention; when work is higher, dGpet should
%become less favorable

HnaughtDA = 10;
%electronic coupling at Van der Waals distance
rnaughtDA = 1.5;
%Van der Waals distance (angstroms)
beta1 = 0.17;
beta2 = 0.215;
%coupling efficiency of the wire in question (angstroms^-1)
HDA1 = HnaughtDA*exp(-beta1*(10*r-rnaughtDA));
HDA2 = HnaughtDA*exp(-beta2*(10*r-rnaughtDA));
%calculates donor-acceptor coupling of electron transfer (with correction on r to be in angstroms); still need to
%poke in literature for best approaches for Van der Waals radius and
%coupling (currently rnaughtDA set to nitrogen radius and HnaughtDA set to
%10 as placeholder (I seem to see this in the 10-100 range). 
%Beta set to minimum value for polystyrene wires as a placeholder

hbar = 6.582119569e-16;
%planck's constant over 2 pi because physicists are nerds; also let's keep
%everything in eV since that is useful
kb = 8.617333262e-5;
%Boltzmann constant in eV
T = 310.15;
%Temperature in Kelvin (assuming T=37 degrees for cells)
lambda = 1;
%solvent reorganization energy; still thinking of the best way to get at
%this (literature? approximations? fitting to known data? currently 1 is a placeholder)
%I seem to see this range from 0.5 to 2ish
kpet1 = sqrt(pi/(hbar*lambda*kb*T)).*(HDA1^2).*exp(-((lambda+dGpet).^2)./(4*lambda*kb*T));
kpet2 = sqrt(pi/(hbar*lambda*kb*T)).*(HDA2^2).*exp(-((lambda+dGpet).^2)./(4*lambda*kb*T));
%at last we have our full kpet calculation from our measured and calculated
%parameters; kpet should be in s^-1

tauprot = 3.5e-9;
%fluorescence lifetime of fully protonated voltagefluor; should be equal to
%1/(kfl+knr) and might be the most tractable experimental way to get at
%those values (currently set to VF2.0Cl lifetime as a placeholder)
kprot = 1/tauprot;
%for variable simplicity later
ref = find(Vmem==-60);
dFoverF_VF = (kprot+kpet1(ref))./(kprot+kpet1)-1;
dFoverF_FVF = (kprot+kpet2(ref))./(kprot+kpet2)-1;

%uses kpet to calculate dF/F for that particular voltage step from the
%desired reference voltage (defaults to -60 mV)


Vmempatch = Vmem(2001:4001);
dFoverFpatch_VF = dFoverF_VF(2001:4001);
dFoverFpatch_FVF = dFoverF_FVF(2001:4001);
Vmempatchint0 = Vmempatch + 60;
%selects the range of Vmem and dF/F actually used in our patch experiments
%(-100 mV to +100 mV). Also shifts dF/F values to set -60 mV to 0 dF/F
dFslope_VF = transpose(Vmempatchint0)\transpose(dFoverFpatch_VF);
dFslope_FVF = transpose(Vmempatchint0)\transpose(dFoverFpatch_FVF);
%finds slope of voltage sensitivity (will be in absolute numbers per 1 mV)
sensitivity_VF = 100*100*dFslope_VF;
sensitivity_FVF = 100*100*dFslope_FVF;
%multiply by 100 to get per 100 mV and then another 100 to get percentage;
%this is the dF/F per 100 mV percent value that we typically report
ratio=sensitivity_VF/sensitivity_FVF;

fprintf('dF/F for VF is %g percent per 100mV \n', sensitivity_VF);
fprintf('dF/F for FVF is %g percent per 100mV \n', sensitivity_FVF);
fprintf('ratio is %g \n', ratio);

dFlinfit_VF = dFslope_VF.*(Vmem+60);
dFlinfit_FVF = dFslope_FVF.*(Vmem+60);
%takes linear fit from above and generates line for plotting

figure(1);
hold on;
plot(Vmem,dFoverF_VF,Vmem,dFoverF_FVF);
plot(Vmem,dFlinfit_VF,Vmem,dFlinfit_FVF);

figure(2);
hold on;
plot(Vmem(2001:200:4001),dFoverF_VF(2001:200:4001),'m.','markersize',25);
plot(Vmem(2001:200:4001),dFoverF_FVF(2001:200:4001),'g.','markersize',25);
plot(Vmem(2001:4001),dFlinfit_VF(2001:4001),'m-');
plot(Vmem(2001:4001),dFlinfit_FVF(2001:4001),'g-');
xticks([-100:20:100]);

figure(3);
hold on;
plot(Vmem(2001:4001),log10(kpet1(2001:4001)));
plot(Vmem(2001:4001),log10(kpet2(2001:4001)));
%plots overall fit as well as linear model from patching range
%plotting can be modified as needed for generation of prettier figures




