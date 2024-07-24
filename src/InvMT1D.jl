# Dieno Diba 2022
# Regularized 1D magnetotelluric inversion
# Gauss-Newton minimization, Laplacian roughening matrix

# Number of layer and thickness of layer are fixed
# The inversion only optimizes the resistivity of each layer 

# ============ INPUT FROM USER ============
# Filename of input data (e.g. "syndata.dat")
finp = "syndata.dat"
# Filename of setting (e.g. "setting.dat")
fset = "setting.dat"
# Filenam e of the output file (e.g. "run01.out" )
fout = "run01.out"
# fout = string(Dates.format(now(),"yyyymmdd_HHMMSS"),".out")   # Using current time 
# ============ INPUT FROM USER ============

# ============== MAIN PART =========== begin
using DelimitedFiles
using Statistics
using SparseArrays
using LinearAlgebra
using Printf
using Dates

# Read inversion setting 
set = readdlm(fset;skipstart=1)
nlayer = set[1,1]
hei = set[2,1]
rst = set[3,1]
λ = set[4,1]
ε = set[5,1]
itmax = set[6,1]

# Read the (synthetic) data
dat = readdlm(finp;skipstart=2)
period = dat[:,1]
nper = length(period)
orh = log10.(dat[:,2])
erh = dat[:,3]./dat[:,2]*log(10)
oph = pi/180*dat[:,4]
eph = pi/180*dat[:,5]
dobs = [orh;oph]
N = nper*2

# Weighting matrix
W = zeros(N,N)
for idN = 1:N
    if idN <= nper
        W[idN,idN] = erh[idN]^-2
    else
        W[idN,idN] = eph[idN-nper]^-2
    end
end
W = W/mean(W)

rho = ones(nlayer)*rst
thi = ones(nlayer-1)*hei
M = nlayer
m0 = log10.(rho)

# Simple Laplacian roughening matrix 
L = zeros(M-2,M)
for idL = 1:M-2
    L[idL,idL+1] = 2
    L[idL,idL] = -1
    L[idL,idL+2] = -1
end

# Forward modeling the starting model
crh,cph = FwdMT1D(rho,thi,1 ./period)
dcal = [log10.(crh);pi/180*cph]
# Objective function of the starting model
phi = (dobs-dcal)'*sparse(W)*(dobs-dcal) .+ λ*m0'*sparse(L)'*sparse(L)*m0;

m = m0
for itr = 2:itmax
    # Compute Jacobian matrix using Delta method
    A = zeros(N,M)
    for idM = 1:M
        rtmp = 10 .^m
        rper = 0.0001
        rtmp[idM] = rtmp[idM] + rper
        prh,pph = FwdMT1D(rtmp,thi,1 ./period)
        dper = [log10.(prh);pi/180*pph]
        A[:,idM] = (dper-dcal)/(log10(rtmp[idM]) - log10(10^m[idM]))
    end
    
    # Update the model using Gauss-Newton method
    if itr < itmax
        grd = -2*(A'*sparse(W)*(dobs-dcal)) .+ 2*λ*sparse(L)'*sparse(L)*m
        hes = 2*(A'*sparse(W)*A) .+ 2*λ*sparse(L)'*sparse(L)
        dmph = hes .+ ε*I(M)
        global m = m - dmph\grd
    end

    # Forward modeling the new model
    frho = 10 .^m
    global crh,cph = FwdMT1D(frho,thi,1 ./period)
    global dcal = [log10.(crh);pi/180*cph]
    # Objective function of the new model
    global phi = vcat(phi,(dobs-dcal)'*sparse(W)*(dobs-dcal) .+ λ*m0'*sparse(L)'*sparse(L)*m0);

    # Terminate the inversion conditionally
    if phi[itr] > phi[itr-1]
        println("Obj. function diverged after ",string(itr)," iterations")
        break
    elseif phi[itr-1] - phi[itr] < 0.01*phi[itr-1]      # stopping crit: 1%
        println("Obj. function converged after ",string(itr)," iterations")
        break
    end
end
# ============== MAIN PART ============ end

# Export the output to a file automatically
fid = open(fout, "w")
print(fid,Dates.format(now(),"yyyymmdd-HHMMSS"),"\n\n")
print(fid,"Calculated vs Observed Data","\n")
print(fid,"   Period(s)  orh(Ohm.m)  crh(Ohm.m)    oph(Deg)    cph(Deg)\n")
for idf = 1:nper
    @printf(fid,"%12.2e%12.2e%12.2e%12.2e%12.2e\n",period[idf],10 .^orh[idf],crh[idf],180/pi*oph[idf],cph[idf])    
end
print(fid,"\n")
print(fid,"Final resistivity structure","\n")
print(fid," No.  Rho(Ohm.m)        h(m)\n")
for idn = 1:nlayer
    if idn < nlayer    
        @printf(fid,"%4i%12.2e%12.2e\n",idn,10^m[idn],thi[idn])
    else
        @printf(fid,"%4i%12.2e\n",idn,10^m[idn])
    end
end
print(fid,"\n")
print(fid,"History of objective function","\n")
print(fid,"   Epoch  Obj.function\n")
for itr = 1:length(phi)
    @printf(fid,"%8i%14.2e\n",itr-1,phi[itr])
end
close(fid)
