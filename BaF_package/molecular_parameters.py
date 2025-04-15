#Constants
uB = 1.399624604;
uN = uB/1836
cmIn2MHz = 0.0299792e6


#138

#Ground state parameters



Y01=6491.3962#6491.3946
Y11=-34.8831#-34.8784
Y21=15.93e-3#13.0288e-3
Y02=-5.5250e-3#-5.5248e-3
Y12=-9.43e-6#-9.7632e-6

gamma_00=80.984
gamma_10=-58.4e-3
delta_gamma = 112e-6
BN=Y01+Y11*1/2+Y21*1/4
DN=-Y02-Y12*1/2
gamma=gamma_00+gamma_10*1/2

bBa = 2303.4; cBa = 75.1965; bFBa = bBa + cBa/3;
eq0Q = -143.6812;
bF = 63.41446; cF = 7.30504; bFF = bF + cF/3;
cI = 0



#Ground state Zeeman parameters
gS = 2.00197
gI2 = 5.258
gI1 = 0.94/1.5
grot = -0.048
gl = -0.00594


#Excited state
Gamma = 2.7 #In MHz. The program should convert it to angualr frequency as necessary
A=632.28165*cmIn2MHz
AD = 0.0310*1e-3*cmIn2MHz
p2q=-0.2578*cmIn2MHz

a_ex = 26.55 #this could be made larger to simulate larger separation betwee th
b_ex = -0.2303
c_ex = -5.3094

h_Ba_12  = 206.6
d_Ba  = 254.3
eq0Q1 = -88.3

h_F_12 = (a_ex-1/2*(b_ex+c_ex))*1
d_F = 3.58*1



Bex=0.21189575*cmIn2MHz
T00=11946.3152*cmIn2MHz

glp = -0.536 ; gLp = 0.98;
