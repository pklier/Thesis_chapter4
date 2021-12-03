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
weight1 = [0.0036    0.0101    0.0165    0.0232    0.0299    0.0362    0.0417    0.0442    0.0480    0.0524    0.0544 0.0553    0.0554    0.0535    0.0522    0.0499    0.0491    0.0466    0.0423    0.0371    0.0320    0.0282 0.0239    0.0203    0.0180    0.0162    0.0138    0.0118    0.0097    0.0074    0.0058    0.0042    0.0033 0.0021    0.0011    0.0006         0];
weight2 = [0.0128800250000000,0.0383919600000000,0.0587751260000000,0.0761683420000000,0.0905182160000000,0.0993812810000000,0.100433417000000,0.0982600500000000,0.0873963570000000,0.0744943470000000,0.0621199750000000,0.0509359300000000,0.0406218590000000,0.0312091710000000,0.0238724870000000,0.0174152010000000,0.0117116830000000,0.00777010100000000,0.00496545200000000,0.00377198500000000,0.00295854300000000,0.00255653300000000,0.00178077900000000,0.000977000000000000,0.000465000000000000,0.000126000000000000,3.45000000000000e-05,9.42000000000000e-06,0,0,0,0,0,0,0,0,0];
%angle between wire and electric field; should be modeled computationally
%or this program can be modified to calculate this from measured dF/F.
%Currently set to the computationally derived distribution of angles for
%VF2.1Cl from Rishi's paper
dmem = 4;
%thickness of the membrane in question (nm); 4 nm is a ballpark number. Not
%dubelco's modified eagle media
w1 = (r.*Vmem./(1000*dmem)).*sum(weight1.*cosd(theta));
w2 = (r.*Vmem./(1000*dmem)).*sum(weight2.*cosd(theta));
%w = (r.*Vmem./(1000*dmem)).*cosd(theta);

%calculates work to move electron in PeT (signs and charge accounted for, should be net positive so watch out
%for cosine range) based on above parameters (takes cosine in degrees);
%1000 in denominator accounts for mV to V conversion
%sanity check: work is positive when Vmem is positive (takes energy to move
%an electron from positive to negative) and work is negative when Vmem is
%negative (electron wants to go from negative to positive)

dGpet1 = dEox - dEred - dEoo + w1;
dGpet2 = dEox - dEred - dEoo + w2;
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
beta = 0.01;
%coupling efficiency of the wire in question (angstroms^-1)
HDA = HnaughtDA*exp(-beta*(10*r-rnaughtDA));
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
kpet1 = sqrt(pi/(hbar*lambda*kb*T)).*(HDA^2).*exp(-((lambda+dGpet1).^2)./(4*lambda*kb*T));
kpet2 = sqrt(pi/(hbar*lambda*kb*T)).*(HDA^2).*exp(-((lambda+dGpet2).^2)./(4*lambda*kb*T));
%at last we have our full kpet calculation from our measured and calculated
%parameters; kpet should be in s^-1

tauprot = 3.5e-9;
%fluorescence lifetime of fully protonated voltagefluor; should be equal to
%1/(kfl+knr) and might be the most tractable experimental way to get at
%those values (currently set to VF2.0Cl lifetime as a placeholder)
kprot = 1/tauprot;
%for variable simplicity later
ref = find(Vmem==-60);
dFoverFsulf = (kprot+kpet1(ref))./(kprot+kpet1)-1;
dFoverFdisulf = (kprot+kpet2(ref))./(kprot+kpet2)-1;
%uses kpet to calculate dF/F for that particular voltage step from the
%desired reference voltage (defaults to -60 mV)


Vmempatch = Vmem(2001:4001);
dFoverFpatchsulf = dFoverFsulf(2001:4001);
dFoverFpatchdisulf = dFoverFdisulf(2001:4001);
Vmempatchint0 = Vmempatch + 60;
%selects the range of Vmem and dF/F actually used in our patch experiments
%(-100 mV to +100 mV). Also shifts dF/F values to set -60 mV to 0 dF/F
dFslopesulf = transpose(Vmempatchint0)\transpose(dFoverFpatchsulf);
dFslopedisulf = transpose(Vmempatchint0)\transpose(dFoverFpatchdisulf);
%finds slope of voltage sensitivity (will be in absolute numbers per 1 mV)
sensitivitysulf = 100*100*dFslopesulf;
sensitivitydisulf = 100*100*dFslopedisulf;
%multiply by 100 to get per 100 mV and then another 100 to get percentage;
%this is the dF/F per 100 mV percent value that we typically report

fprintf('dF/F of monosulf is %g percent per 100mV \n', sensitivitysulf);
fprintf('dF/F of disulf is %g percent per 100mV \n', sensitivitydisulf);

dFlinfitsulf = dFslopesulf.*(Vmem+60);
dFlinfitdisulf = dFslopedisulf.*(Vmem+60);
%takes linear fit from above and generates line for plotting


%figure(1);
%hold on;
%plot(Vmem,dFoverFsulf);
%plot(Vmem,dFlinfit);

figure(2);
hold on;
plot(Vmem(2001:200:4001),dFoverFsulf(2001:200:4001),'m.','markersize',25);
plot(Vmem(2001:4001),dFlinfitsulf(2001:4001),'m-');

plot(Vmem(2001:200:4001),dFoverFdisulf(2001:200:4001),'g.','markersize',25);
plot(Vmem(2001:4001),dFlinfitdisulf(2001:4001),'g-');

xticks([-100:20:100]);

%figure(3);
%hold on;
%plot(Vmem(2001:4001),log10(kpet1(2001:4001)));
%plots overall fit as well as linear model from patching range
%plotting can be modified as needed for generation of prettier figures




