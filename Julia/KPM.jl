__precompile__()

module KPM

using LinearAlgebra, FFTW, LinearMaps, Arpack

export DOS_KPM, computeDOS, ρ0, computeDOSslow, normalizeH, normalizeHsymmetric, normalizeH_coeffs

function DensityOfStates!(Evals::Vector{Vector{Float64}},ρEtemp::Vector{Vector{Float64}},
    ρE::Vector{Vector{Float64}},ρE²::Vector{Vector{Float64}},α0::Vector{ComplexF64},
    α1::Vector{ComplexF64},α2::Vector{ComplexF64},μ::Vector{Float64},H,
    NR::Integer,intermediateNCs::Vector{Int64},intermediateNtildes::Vector{Int64},
    intermediate_ρ0s::Vector{Float64},intermediate_ρ0²s::Vector{Float64},
    intermediate_d2ρ0s::Vector{Float64},intermediate_d2ρ0²s::Vector{Float64},
    intermediate_d4ρ0s::Vector{Float64},intermediate_d4ρ0²s::Vector{Float64};
    Emin::Float64=-1.0,Emax::Float64=1.0,Kernel::Function=JacksonKernel,
    intermediateNCs0::Vector{Int64}=Array{Int64}(undef,0), outputBW::Bool=false)

    dimrange, dim = size(H)
    # @assert dimrange==dim "H must be square"
    # @assert length(ρEtemp)==length(Evals) && length(ρEtemp)==length(ρE) && length(ρEtemp)==length(ρE²) "all E and ρ(E) must be same size"
    # @assert length(α0)==dim && length(α1)==dim && length(α2)==dim "α vectors must be dimension of H"
    NC = length(μ)
    Ntilde = length(ρEtemp)

    a,b = normalizeH_coeffs(H)
    if outputBW
        println("$(a) $(b)")
    end
    # μ .= 0.0
    computeμs_randomphase!(μ,α0,α1,α2,H,a,b,div(NC,2),NR,dim)


    for i=eachindex(intermediateNCs)
        iNC = intermediateNCs[i]
        iNtilde = intermediateNtildes[i]
        computeDOSslow!(Evals[i],ρEtemp[i],μ[1:iNC],a,b,iNC,iNtilde;Emin=Emin,Emax=Emax,Kernel=Kernel)

        for j=eachindex(ρEtemp[i])
            @inbounds ρE[i][j] += ρEtemp[i][j]
            @inbounds ρE²[i][j] += ρEtemp[i][j]^2
        end

        rho0, d2rho0, d4rho0 = ρ0(μ[1:iNC], a, b)
        intermediate_ρ0s[i] += rho0
        intermediate_d2ρ0s[i] += d2rho0
        intermediate_d4ρ0s[i] += d4rho0
        intermediate_ρ0²s[i] += rho0^2
        intermediate_d2ρ0²s[i] += d2rho0^2
        intermediate_d4ρ0²s[i] += d4rho0^2

    end

    maxNCind = length(intermediateNCs)

    for i=eachindex(intermediateNCs0)
        iNC = intermediateNCs0[i]
        rho0, d2rho0, d4rho0 = ρ0(μ[1:iNC], a, b)
        intermediate_ρ0s[maxNCind+i] += rho0
        intermediate_d2ρ0s[maxNCind+i] += d2rho0
        intermediate_d4ρ0s[maxNCind+i] += d4rho0
        intermediate_ρ0²s[maxNCind+i] += rho0^2
        intermediate_d2ρ0²s[maxNCind+i] += d2rho0^2
        intermediate_d4ρ0²s[maxNCind+i] += d4rho0^2
    end

end

function DensityOfStates!(Evals::Vector{Float64},ρEtemp::Vector{Float64},
    ρE::Vector{Float64},ρE²::Vector{Float64},α0::Vector{ComplexF64},
    α1::Vector{ComplexF64},α2::Vector{ComplexF64},μ::Vector{Float64},H,
    NR::Integer;Emin::Float64=-1.0,Emax::Float64=1.0,Kernel::Function=JacksonKernel,
    outputBW::Bool=false)

    dimrange, dim = size(H)
    # @assert dimrange==dim "H must be square"
    # @assert length(ρEtemp)==length(Evals) && length(ρEtemp)==length(ρE) && length(ρEtemp)==length(ρE²) "all E and ρ(E) must be same size"
    # @assert length(α0)==dim && length(α1)==dim && length(α2)==dim "α vectors must be dimension of H"
    NC = length(μ)
    Ntilde = length(ρEtemp)

    a,b = normalizeH_coeffs(H)
    μ .= 0.0
    computeμs_randomphase!(μ,α0,α1,α2,H,a,b,div(NC,2),NR,dim)
    computeDOSslow!(Evals,ρEtemp,μ,a,b,NC,Ntilde;Emin=Emin,Emax=Emax,Kernel=Kernel)

    for i=eachindex(ρEtemp)
        @inbounds ρE[i] += ρEtemp[i]
        @inbounds ρE²[i] += ρEtemp[i]^2
    end

    if outputBW
        return ρ0(μ,a,b), a, b
    else
        return ρ0(μ,a,b)
    end
end

function DensityOfStates(H,NC::Integer,NR::Integer,Ntilde::Integer; outputBW::Bool=false,
    Emin::Float64=0.0, Emax::Float64=0.0, Kernel::Function=JacksonKernel)

    dimrange, dim = size(H)
    @assert dimrange==dim "H must be square"
    a, b = normalizeH_coeffs(H)
    μ::Vector{Float64} = zeros(Float64,NC)
    α0 = Array{ComplexF64}(undef,dim)
    α1 = Array{ComplexF64}(undef,dim)
    α2 = Array{ComplexF64}(undef,dim)

    Evals = Array{Float64}(undef,Ntilde)
    rhoE = Array{Float64}(undef,Ntilde)

    computeμs_randomphase!(μ,α0,α1,α2,H,a,b,div(NC,2),NR,dim)
    if Emin!=Emax
        computeDOSslow!(Evals,rhoE,μ,a,b,NC,Ntilde;Emin=Emin,Emax=Emax,Kernel=Kernel)
    else
        computeDOS!(rhoE,Evals,μ,a,b,NC,Ntilde; Kernel=Kernel)
    end

    if outputBW
        return Evals, rhoE, a, b
    else
        return Evals, rhoE
    end
    # return Evals, rhoE

end

DOS(H,a::Float64,b::Float64,NC::Integer,NR::Integer) = DOS_KPM(H,a,b,div(NC,2),NR,size(H)[1])
DOS(H,NC::Integer,NR::Integer) = DOS_KPM(H,div(NC,2),NR,size(H)[1])
DOS_KPM(H,NC::Integer,NR::Integer) = DOS_KPM(H,NC::Integer,NR::Integer, size(H)[1])
DOS_KPM(H,NC::Integer,NR::Integer,dim::Integer) = DOS_KPM(H,1.0,0.0,NC,NR,dim)

function DOS_KPM(H,a::Float64,b::Float64,NC::Integer,NR::Integer,dim::Integer)::Vector{Float64}

    dimrange, dim = size(H)
    @assert dimrange == dim "H must be square"

    μ::Array{Float64} = zeros(Float64,2NC)
    α0 = Array{ComplexF64}(undef,dim)
    α1 = Array{ComplexF64}(undef,dim)
    α2 = Array{ComplexF64}(undef,dim)

    computeμs_randomphase!(μ,α0,α1,α2,H,a,b,NC,NR,dim)

    return μ
end

function computeμs_randomphase!(μ::Vector{Float64},α0::Vector{ComplexF64},
    α1::Vector{ComplexF64},α2::Vector{ComplexF64},H,a::Float64,b::Float64,
    NC::Integer,NR::Integer,dim::Integer)

    μ .= 0.0

    for r=1:NR
        for j=eachindex(α0)
            @inbounds α0[j] = exp(2im*pi*rand(Float64))/sqrt(dim)
        end
        symmetricKPM!(μ,α0,α1,α2,H,a,b,NC)
    end
    μ ./= NR

    maximum(abs.(μ))<=1.0 || error("bounds computed incorrectly")

end

function symmetricKPM!(μ::Vector{Float64},α0::Vector{ComplexF64},
    α1::Vector{ComplexF64},α2::Vector{ComplexF64},H,a::Float64,b::Float64,
    NC::Integer)

    mu1::Float64 = 0.0

    mul!(α1,H,α0)
    @. α1 = (α1 - b*α0)/a

    μ[1] += 1.0
    mu1 = real(dot(α0,α1))
    μ[2] += mu1

    for n=2:NC

        mul!(α2,H,α1)
        @. α2 = 2*(α2-b*α1)/a - α0

        @inbounds μ[2n-1] += 2*real(dot(α1,α1)) - 1.0
        @inbounds μ[2n] += 2*real(dot(α2,α1)) - mu1

        α2, α0, α1 = (α0, α1, α2)

    end

end


JacksonKernel(n::Integer,N::Integer) = 1/(N+1)*((N+1-n)*cos(pi*n/(N+1))
            + sin(pi*n/(N+1))*cot(pi/(N+1)))

LorentzKernel(n,N::Integer,λ::Float64) = sinh(λ*(1-n/N))/sinh(λ)
FejérKernel(n,N::Integer) = 1 - n/N
LanczosKernel(n,N,M::Integer) = (sin(π*n/N)/(π*n/N))^M
WangZungerKernel(n,N,α::Float64,β::Float64) = exp(-(α*n/N)^β)
DirichletKernel(n,N) = 1.0

function computeDOS(μ::Array{Float64,1},a::Float64,b::Float64,NC::Integer,
    Ntilde::Integer; Kernel::Function=JacksonKernel)

    @assert Ntilde>=NC "Ntilde must be >= NC"
    μtilde = zeros(Float64,Ntilde)
    xs = Array{Float64}(undef,Ntilde)
    computeDOS!(μtilde,xs,μ,a,b,NC,Ntilde;Kernel=Kernel)
    return xs, μtilde
end

    #computeDOS!(rhoE,Evals,μ,a,b,NC,Ntilde; Kernel=Kernel)
function computeDOS!(μtilde::Vector{Float64},xs::Vector{Float64},
    μ::Vector{Float64},a::Float64,b::Float64,NC::Integer,Ntilde::Integer;
    Kernel::Function=JacksonKernel)

    @. μtilde = 0.0
    for i=1:NC
        @inbounds μtilde[i] = μ[i]*Kernel(i-1,NC)
    end
    μtilde[1] = μtilde[1]/sqrt(2)
    for j=eachindex(xs)
        @inbounds xs[j] = cos(pi*(j-0.5)/Ntilde)
    end
    idct!(μtilde) # Needs memory allocations
    @. μtilde *= sqrt(2*Ntilde)/(a*pi*sqrt(1-xs^2))
    @. xs = a*xs + b
end

# chebyshevT(n::Integer,x::Float64)::Float64 = chebyshevT(n,x)
chebyshevT(n::Integer,x::Float64)::Float64 = cos(n*acos(x))
chebyshevU(n,x) = sin((n+1)*acos(x))/sqrt(1 - x^2)

function computeDOSslow(μ::Array{Float64,1},a::Float64,b::Float64,NC::Integer,
    Ntilde::Integer;Emin::Float64=-1.0,Emax::Float64=1.0,
    Kernel::Function=JacksonKernel)

    Evals = Array{Float64}(undef,Ntilde)
    rhoE = Array{Float64}(undef,Ntilde)

    computeDOSslow!(Evals,rhoE,μ,a,b,NC,Ntilde;Emin=Emin,Emax=Emax,Kernel=Kernel)

    return Evals, rhoE
end


#    computeDOSslow!(Evals,rhoE,μ,a,b,NC,Ntilde;Emin=Emin,Emax=Emax,Kernel=Kernel)
function computeDOSslow!(Evals::Vector{Float64},rhoE::Vector{Float64},
    μ::Array{Float64,1},a::Float64,b::Float64,NC::Integer,
    Ntilde::Integer;Emin::Float64=-1.0,Emax::Float64=1.0,
    Kernel::Function=JacksonKernel)

    # for i=1:Ntilde
    for i=eachindex(Evals)
        @inbounds Evals[i]= Emin+(Emax-Emin)*(i-1)/(Ntilde-1)
    end

    for k=eachindex(rhoE)
        @inbounds rhoE[k] = μ[1]*Kernel(0,NC)
        for ii=2:NC
            @inbounds rhoE[k] += 2*μ[ii]*(Kernel(ii-1,NC)*
              chebyshevT(ii-1,(Evals[k]-b)/a))
        end
    end
    @. rhoE /= a*pi*sqrt( 1 - (Evals - b)^2/a^2)

end

function ρ0Symmetric(μ::Array{Float64,1},a::Float64)
    NC::Integer = length(μ)
    ρ₀::Float64 = μ[1]*JacksonKernel(0,NC)/pi
    ∂²ρ₀ = ρ₀
    ∂⁴ρ₀ = 9ρ₀
    for n=1:(div(NC,2)-1)
        newterm = 2/pi*μ[2n+1]*JacksonKernel(2n,NC)*(-1)^n
        ρ₀ += newterm
        ∂²ρ₀ += newterm*(1 - (2n)^2)
        ∂⁴ρ₀ += newterm*(9 - 10*(2n)^2+(2n)^4)
    end
    return ρ₀/a, ∂²ρ₀/a^3, ∂⁴ρ₀/a^5
end

function ρ0(μ::Array{Float64,1},a::Float64,b::Float64; E0::Float64=0.0)
    NC = length(μ)
    x0 = (E0-b)/a
    ρ₀::Float64 = μ[1]*JacksonKernel(0,NC)/(pi*sqrt(1 - x0^2))
    ∂²ρ₀::Float64 = ρ₀*(1 + 2x0^2)/(1 - x0^2)^2
    ∂⁴ρ₀::Float64 = ρ₀*(9 + 24*x0^2*(3+ x0^2))/(1 - x0^2)^(4)
    for n=1:(NC-1)
        @inbounds newterm = 2/pi*μ[n+1]*JacksonKernel(n,NC)/sqrt(1 - x0^2)
        ρ₀ += newterm*chebyshevT(n,x0)
        ∂²ρ₀ += newterm*( (1 + 2x0^2 - n^2(1 - x0^2))/(1 - x0^2)*chebyshevT(n,x0)
                           + 3*n*x0/(1 - x0^2)*chebyshevU(n-1,x0))
        ∂⁴ρ₀ += newterm*( 3*(3+ x0^2*(n^2*(n^2 - 25)*(1-x0^2) + 8*(3+x0^2)))/(1 - x0^2)^4*chebyshevT(n,x0)
                         + n*x0*(55 - 39n - 10n^2 + 3n^3 - 2*(n-5)*(5+5n+2n^2)*x0^2)/(1-x0^2)^3*chebyshevU(n-1,x0)
                         + n^2*(n^2-10)/(1 - x0^2)^3*chebyshevU(n-4,x0) )
    end
    return ρ₀/a, ∂²ρ₀/a^3, ∂⁴ρ₀/a^5
end

function normalizeHsymmetric(H;ϵ::Float64=0.01)
    es, vs = eigs(H)
    Emax = maximum(abs.(es))
    Emin = -Emax
    a = (Emax - Emin)/(2 - ϵ)
    H = H/a
    return a, H
end

function normalizeH(H;ϵ::Float64=0.01,maxiter::Int64=1000,tol::Float64=1e-2)
    es, vs = eigs(H; maxiter = maxiter,tol = tol)
    Emax = maximum(real.(es))
    Emin = minimum(real.(es))
    if Emin >= 0.0
        es, vs = eigs(H - Emax*I; maxiter = maxiter,tol = tol)
        Emin = minimum(real.(es)) + Emax
    elseif Emax<=0.0
        es, vs = eigs(H - Emin*I; maxiter = maxiter,tol = tol)
        Emax = maximum(real.(es)) + Emin
    end
    a = (Emax - Emin)/(2 - ϵ)
    b = (Emax + Emin)/2
    H = (H-b*I)/a
    return a, b, H
end

function normalizeH_coeffs(H;ϵ::Float64=0.01,maxiter::Int64=1000,tol::Float64=1e-2)
    es = real.(getindex(eigs(H; maxiter = maxiter,tol = tol),1))
    Emax = maximum(es)
    Emin = minimum(es)
    if Emin >= 0.0
        es = real.(getindex(eigs(H - Emax*I; maxiter = maxiter,tol = tol),1))
        Emin = minimum(es) + Emax
    elseif Emax<=0.0
        es, vs = real.(getindex(eigs(H - Emin*I; maxiter = maxiter,tol = tol),1))
        Emax = maximum(es) + Emin
    end
    a = (Emax - Emin)/(2 - ϵ)
    b = (Emax + Emin)/2

    return a, b
end

end
