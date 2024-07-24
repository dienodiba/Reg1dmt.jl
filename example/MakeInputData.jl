# Dieno Diba 2022
# Make synthetic data for the input of 1D MT inversion (REG1DMT.jl)

# ============ INPUT FROM USER ============
# Resistivity (Ohm.m)
rho = [100;1;100]
# Layer thickness (m)
h = [50,40]
# Frequency
fmin = -3   # min in logscale
fmax = 4    # max in logscale
fpts = 20   # number of points
# Output Filename (e.g. "syndata.dat")
fout = "syndata.dat"
# ============ INPUT FROM USER ============

# Forward Computation
freq = 10 .^ range(fmax, stop=fmin, length=fpts);
nfreq = length(freq)
period = 1 ./ freq

rha, pha = FwdMT1D(rho,h,freq)

# Add some gaussian noise to the synthetic data
srh = randn(nfreq);
drh = 0.1*srh;
rhn = 10 .^(log10.(rha) + drh);
erh = abs.(drh).*rhn/log(10);
dph = 3*srh;
phn = pha + dph;
eph = abs.(dph);

# Export to a file automatically
using Printf
fid = open(fout, "w")
print(fid,"Nperiod = ",nfreq,"\n")
print(fid,"   Period(s)  rha(Ohm.m)  erh(Ohm.m)    pha(deg)    eph(deg)\n")
for idf = 1:nfreq
    @printf(fid,"%12.2e%12.2e%12.2e%12.2e%12.2e\n",period[idf],rhn[idf],erh[idf],phn[idf],eph[idf])    
end
close(fid)