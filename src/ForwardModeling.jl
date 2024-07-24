# Dieno Diba 2022
# Magnetotelluric 1D forward modeling

function ForwardModeling(rho,h,freq)
    mu = 4*pi*1e-7
    omega = 2*pi*freq
    nfreq = length(freq)
    nrho = length(rho)
    rha = zeros(nfreq,1);
    pha = zeros(nfreq,1);

    for idf = 1:nfreq
        Z = zeros(ComplexF64,nrho)
        Z0 = zeros(ComplexF64,nrho-1); 
        R = zeros(ComplexF64,nrho-1); 
        k = zeros(ComplexF64,nrho-1);
        Z[nrho] = sqrt(im*omega[idf]*mu*rho[nrho])
        for idl = nrho-1:-1:1
            Z0[idl] = sqrt(im*omega[idf]*mu*rho[idl])
            R[idl] = (Z0[idl] - Z[idl+1]) / (Z0[idl] + Z[idl+1])
            k[idl] = Z0[idl] / rho[idl]
            Z[idl] = Z0[idl] * ((1 - (R[idl]*exp(-2*k[idl]*h[idl]))) / (1 + (R[idl]*exp(-2*k[idl]*h[idl]))))
        end
        rha[idf] = (1/(omega[idf]*mu)) * (abs(Z[1]))^2
        pha[idf] = 180/pi*angle(Z[1])
    end
    
    return rha, pha
end
