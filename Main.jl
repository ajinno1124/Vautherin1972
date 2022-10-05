include("SkyrmeParams.jl")
include("MyLib.jl")
using Parameters
using LinearAlgebra
using .NuclParameters
using .LambdaParameters
using .MyLib

@consts begin
    # Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
    mpMeV=938.27208816
    mnMeV=939.5654205
    mAveMeV=(mpMeV+mnMeV)/2
    mΛMeV=1115.683 #\pm 0.006 MeV
    ħc=197.3269804
    e2MeVfm=1.4400

    Nmesh=600
    Nmatch=150
    rmax=30
    lmax=7
end

mutable struct QuantumNumber
    j::Float64
    l::Int64
    B::Int64 #1::p, 2::n, 3::Λ
end

mutable struct SingleParticleState
    QN::QuantumNumber
    E::Float64
    ψ::Vector{Float64}
end

mutable struct AtomNum
    Z::Int64
    N::Int64
    Λ::Int64
end

function getrmesh()
    #Be careful for the SolvePoissionEq
    h=rmax/Nmesh
    rmesh=0:h:h*(Nmesh-1)
    return rmesh
end

function getmass(QN::QuantumNumber)
    mass=0.0
    if QN.B==1
        mass+=mpMeV
    elseif QN.B==2
        mass+=mnMeV
    elseif QN.B==3
        mass+=mΛMeV
    end
    return mass
end

#Initail Condition (Woods-Saxon)
function InitPot(AN::AtomNum,rmesh)
    N=AN.N
    Z=AN.Z
    A=N+Z
    Vcoef=zeros(Float64,3)
    Vcoef[1]=-51-33*((N-Z)/A) #proton
    Vcoef[2]=-52+33*((N-Z)/A) #neutron
    Vcoef[3]=-52*(2/3) #Lambda

    r0=1.25
    a=0.65
    R=r0*A^(1/3)
    h2m=zeros(Float64,(3,Nmesh))
    dh2m=zeros(Float64,(3,Nmesh))
    ddh2m=zeros(Float64,(3,Nmesh))
    V=zeros(Float64,(3,Nmesh))
    W=zeros(Float64,(3,Nmesh))

    for b in 1:3
        if b==1 #proton
            @. h2m[b,:]+=ħc^2/(2*mpMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=(0.22*Vcoef[b]*r0^2/a)*(exp((rmesh[:]-R)/a)/(1+exp((rmesh[:]-R)/a))^2)
        elseif b==2
            @. h2m[b,:]+=ħc^2/(2*mnMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=(0.22*Vcoef[b]*r0^2/a)*(exp((rmesh[:]-R)/a)/(1+exp((rmesh[:]-R)/a))^2)
        elseif b==3
            @. h2m[b,:]+=ħc^2/(2*mΛMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=0
        end
    end

    return h2m,dh2m,ddh2m,V,W
end

#######################################
# Differential Equation AR''+BR+CR=ER #
# R=u, R'=v                           #
#######################################
function CalcABC(QN,h2mB,dh2mB,ddh2m,VB,WB,rmesh)
    j=QN.j
    l=QN.l
    A=zeros(Float64,Nmesh)
    C=zeros(Float64,Nmesh)
    @. A[:] += -h2mB[:]

    C[1]=NaN #singular behavior
    @. C[2:Nmesh] += -0.25*dh2mB[2:Nmesh]^2/h2mB[2:Nmesh]
    @. C[2:Nmesh] += 0.5*ddh2m[2:Nmesh]
    @. C[2:Nmesh] += h2mB[2:Nmesh]*l*(l+1)/(rmesh[2:Nmesh]^2) + VB[2:Nmesh]
    @. C[2:Nmesh] += dh2mB[2:Nmesh]/rmesh[2:Nmesh] + WB[2:Nmesh]/rmesh[2:Nmesh]*(j*(j+1)-l*(l+1)-0.75)

    return A,C
end

#Diff. eq. ψ''(r)+f(r)ψ(r)=0
#Calc ψ[1]=ψ(r-h) by ψ[2]=ψ(r) and ψ[3]=ψ(r+h) using Numerov Method
function Numerov6(ψ::Vector{Float64},f::Vector{Float64},h)
	val=0.0
	val+=(2-5*h^2*f[2]/6)*ψ[2]
	val-=(1+h^2*f[1]/12)*ψ[1]
	val/=(1+h^2*f[3]/12)
	return val
end

function BoundCond(QN,E,A,C,rmesh)
    @assert E<=0
    l=QN.l
    mass=getmass(QN)
    h=rmesh[2]-rmesh[1]
    Rin=zeros(Float64,3)
    Rin[1]=0 # r=rmesh[1]=0
    Rin[2]=h # r=rmesh[2]=h
    if l!=1
        Rin[3]+=(2 - 5*h^2*(C[2]-E)/A[2]/6)*Rin[2]
        #Rin[3]-=0
        Rin[3]/=(1+h^2*(C[3]-E)/A[3]/12)
    else
        Rin[3]+=(2 - 5*h^2*(C[2]-E)/A[2]/6)*Rin[2]
        Rin[3]-=h^2 /12 * (ħc^2/(2*mass)*l*(l+1)*Rin[2]/rmesh[2]^2-E)/A[1]
        Rin[3]/=(1+h^2*(C[3]-E)/A[3]/12)
    end

    Rout=zeros(Float64,3)
    Rout[3]=exp(-(-2*mass/ħc^2*E)^(0.5) * rmesh[Nmesh])
    Rout[2]=exp(-(-2*mass/ħc^2*E)^(0.5) * rmesh[Nmesh-1])
    Rout[1]=exp(-(-2*mass/ħc^2*E)^(0.5) * rmesh[Nmesh-2])

    return Rin,Rout
end



function WronskyEuler(E,QN::QuantumNumber,A,C,rmesh)
    h=rmesh[2]-rmesh[1]

    Rin=zeros(Float64,5)
    Rout=zeros(Float64,5)
    Rin[3:5],Rout[1:3]=BoundCond(QN,E,A,C,rmesh)

    for i in 3:Nmatch+1
        Rin[1:4]=Rin[2:5]
        ψvec=[Rin[3],Rin[4]]
        fvec=[(C[i-1]-E)/A[i-1], (C[i]-E)/A[i], (C[i+1]-E)/A[i+1]]
        Rin[5]=Numerov6(ψvec,fvec,h)
    end

    for i in Nmesh-2:-1:Nmatch-1
        Rout[2:5]=Rout[1:4]
        ψvec=[Rout[3],Rout[2]]
        fvec=[(C[i+1]-E)/A[i+1], (C[i]-E)/A[i], (C[i-1]-E)/A[i-1]]
        Rout[1]=Numerov6(ψvec,fvec,-h)
    end

    dRin=MyLib.diff1st5pt(h,Rin)
    dRout=MyLib.diff1st5pt(h,Rout)

    return Rin[3]*dRout-Rout[3]*dRin
end

##########################################################
#Calculate States by given A, B, C
function NormFact(rmesh,ψ)
    ans=MyLib.IntTrap(rmesh,@. ψ[:]^2)
    ans=sqrt(ans)
    return ans
end

function RadWaveFunc(E,QN::QuantumNumber,A,C,rmesh)
    h=rmesh[2]-rmesh[1]

    R=zeros(Float64,Nmesh)
    R[1:3],R[Nmesh-2:Nmesh]=BoundCond(QN,E,A,C,rmesh)

    for i in 3:Nmatch-1
        ψvec=[R[i-1],R[i]]
        fvec=[(C[i-1]-E)/A[i-1], (C[i]-E)/A[i], (C[i+1]-E)/A[i+1]]
        R[i+1]=Numerov6(ψvec,fvec,h)
    end
    R[1:Nmatch]/=R[Nmatch]

    for i in Nmesh-2:-1:Nmatch+1
        ψvec=[R[i+1],R[i]]
        fvec=[(C[i+1]-E)/A[i+1], (C[i]-E)/A[i], (C[i-1]-E)/A[i-1]]
        R[i-1]=Numerov6(ψvec,fvec,-h)
    end
    R[Nmatch:Nmesh]/=R[Nmatch]

    @. R[:]*=(-A[:])^(-0.5)

    Norm=NormFact(rmesh,R)
    R*=sign(R[2]-R[1])/Norm

    return R
end

function CalcStates(QN::QuantumNumber,h2mB,dh2mB,ddh2mB,VB,WB,rmesh)
    States=SingleParticleState[]

    A,C=CalcABC(QN,h2mB,dh2mB,ddh2mB,VB,WB,rmesh)
    Erange=-100.0:0.0 #In BCS approx., E>0 state also needs to be calculated.
    args=[QN,A,C,rmesh]

    for i in 1:(length(Erange)-1)
        Eans=MyLib.MyBisect(Erange[i],Erange[i+1],WronskyEuler,args,rtol=1e-6) #WronskyEuler(E,QN::QuantumNumber,A,B,C,rmesh)
        if isnan(Eans)==true
            continue
        else
            ψ=RadWaveFunc(Eans,args...)
            push!(States,SingleParticleState(QN,Eans,ψ[:]))
        end
    end

    return States
end

function CalcAllStates(h2m,dh2m,ddh2m,V,W,rmesh)
    AllStates=[SingleParticleState[],SingleParticleState[],SingleParticleState[]]
    for b in 1:3
        States=SingleParticleState[]
        for l in 0:lmax
            for j in max(l-0.5,0.5):l+0.5
                QN=QuantumNumber(j,l,b)
                States=vcat(States,CalcStates(QN,h2m[b,:],dh2m[b,:],ddh2m[b,:],V[b,:],W[b,:],rmesh))
            end
        end
        sort!(States, by=x->x.E)
        AllStates[b]=vcat(AllStates[b],States)
    end

    return AllStates
end

function Calc_occ(AN::AtomNum,AllStates;BCS=false)
    Allocc=[Float64[],Float64[],Float64[]]
    Remain=[AN.Z,AN.N,AN.Λ]
    if BCS==false
        for b in 1:3
            for i=eachindex(AllStates[b])
                j=AllStates[b][i].QN.j
                if Remain[b]>=(2*j+1)
                    push!(Allocc[b],1.0)
                    Remain[b]-=2*j+1
                else
                    push!(Allocc[b],Remain[b]/(2*j+1))
                    Remain[b]=0
                end
            end
        end
    end

    return Allocc
end

function InitialCondition(AN::AtomNum)
    rmesh=getrmesh()
    h2m,dh2m,ddh2m,V,W=InitPot(AN,rmesh)

    InitState=CalcAllStates(h2m,dh2m,ddh2m,V,W,rmesh)
    #Initocc=Calc_occ(AN,InitState)

    return InitState
end

##########################################################
# interpolate the even function value at r=0
function InterPolEvenFunc0(f2,f3,f4)
    val=0.0
    val+=37.0/24.0*f2
    val+=-2.0/3.0*f3
    val+=1.0/8.0*f4

    return val
end

# defene Density, Potential
function Calc_ρ(occ::Vector{Float64},States::Vector{SingleParticleState},rmesh)
    ρ=zeros(Float64,Nmesh)
    for i=eachindex(occ)
        j=States[i].QN.j
        l=States[i].QN.l
        R=States[i].ψ
        if l==0
            R1r1=InterPolEvenFunc0(R[2]/rmesh[2],R[3]/rmesh[3],R[4]/rmesh[4])
            ρ[1]+=occ[i]*(2*j+1)/(4*π)*(R1r1)^2
        end
        @. ρ[2:Nmesh]+=occ[i]*(2*j+1)/(4*π)*(States[i].ψ[2:Nmesh]/rmesh[2:Nmesh])^2
    end
    return ρ
end


function Calc_dρ(ρ,rmesh)
    h=rmesh[2]-rmesh[1]
	P=1 #ρ is even function
    dρ=MyLib.diff1st5pt(h,ρ,P)

    return dρ
end

function Calc_ddρ(ρ,rmesh)
    h=rmesh[2]-rmesh[1]
	P=1 #ρ is even function
    ddρ=MyLib.diff2nd5pt(h,ρ,P)

    return ddρ
end

function Calc_Lapρ(ρ::Vector{Float64},rmesh)
    Lapρ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dρ=Calc_dρ(ρ,rmesh)
    ddρ=Calc_ddρ(ρ,rmesh)
    dρr=zeros(Float64,Nmesh)
    #dρr[1]=dρ[2]/rmesh[2]
    dρr[1]=InterPolEvenFunc0(dρ[2]/rmesh[2],dρ[3]/rmesh[3],dρ[4]/rmesh[4])
    @. dρr[2:Nmesh]=dρ[2:Nmesh]/rmesh[2:Nmesh]
    @. Lapρ[:]+=2*dρr[:] + ddρ[:]
    return Lapρ
end

#Calc_Rr, Calc_dRrを導入しても良いかも
function Calc_τ(occ,States::Vector{SingleParticleState},rmesh)
    τ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dRr=zeros(Float64,Nmesh)
    Rr=zeros(Float64,Nmesh)
    for i=eachindex(occ)
        j=States[i].QN.j
        l=States[i].QN.l
        P_Rr=(-1)^l

        R=States[i].ψ
        if l==0
            Rr[1]=InterPolEvenFunc0(R[2]/rmesh[2],R[3]/rmesh[3],R[4]/rmesh[4])
        else
            Rr[1]=0
        end
        @. Rr[2:Nmesh]=R[2:Nmesh]/rmesh[2:Nmesh]

		dRr=MyLib.diff1st5pt(h,Rr,P_Rr)

        @. τ[:]+=occ[i]*(2*j+1)/(4*π)*dRr[:]^2

        if l==1
            Rr2r2_r0=InterPolEvenFunc0(Rr[2]^2/rmesh[2]^2,Rr[3]^2/rmesh[3]^2,Rr[4]^2/rmesh[4]^2)
            τ[1]+=occ[i]*(2*j+1)/(4*π)*l*(l+1)*Rr2r2_r0
        end
        if l>0
            @. τ[2:Nmesh]+=occ[i]*(2*j+1)/(4*π)*l*(l+1)*Rr[2:Nmesh]^2/rmesh[2:Nmesh]^2
        end
    end
    return τ
end

# τ: Parity even
function Calc_dτ(τ,rmesh)
    h=rmesh[2]-rmesh[1]
    P_τ=1
    dτ=MyLib.diff1st5pt(h,τ,P_τ)

    return dτ
end

function Calc_J(occ,States::Vector{SingleParticleState},rmesh)
    J=zeros(Float64,Nmesh)
    for i=eachindex(occ)
        j=States[i].QN.j
        l=States[i].QN.l
        if l>0
            @. J[2:Nmesh]+=occ[i]*(2*j+1)/(4*π)*(j*(j+1)-l*(l+1)-0.75)*(States[i].ψ[2:Nmesh]/rmesh[2:Nmesh])^2/rmesh[2:Nmesh]
        end
    end
    return J
end

function Calc_dJ(J,rmesh)
    dJ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    P=-1 # J: Parity odd
    dJ=MyLib.diff1st5pt(h,J,P)

    return dJ
end

function Calc_divJ(J::Vector{Float64},rmesh)
    divJ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dJ=Calc_dJ(J,rmesh)
    Jr=zeros(Float64,Nmesh)
    Jr[1]=InterPolEvenFunc0(J[2]/rmesh[2],J[3]/rmesh[3],J[4]/rmesh[4])
    @. Jr[2:Nmesh]=J[2:Nmesh]/rmesh[2:Nmesh]
    @. divJ+=dJ[:]+2*Jr[:]
    return divJ
end

##################################################3
# define potential
function Calc_h2mN(AN::AtomNum,b,aN,aL,ρN::Vector{Float64},ρq::Vector{Float64},ρΛ::Vector{Float64},)
	N=AN.N
	Z=AN.Z
	Λ=AN.Λ
    QN=QuantumNumber(0,0,b)
    m=getmass(QN)
    return @. ħc^2/(2*m)-ħc^2/(2*(mpMeV*Z+mnMeV*N+mΛMeV*Λ))+aN[5]*ρN[:]+aN[6]*ρq[:]+aL[2]*ρΛ[:]
end

function Calc_h2mΛ(AN::AtomNum,aL,ρN::Vector{Float64})
	N=AN.N
	Z=AN.Z
	Λ=AN.Λ
    return @. ħc^2/(2*mΛMeV)-ħc^2/(2*(mpMeV*Z+mnMeV*N+mΛMeV*Λ))+aL[2]*ρN
end

function Calc_VΛΛ(aL,γ,ρN::Vector{Float64},LapρN::Vector{Float64},τN::Vector{Float64},ρp::Vector{Float64},ρn::Vector{Float64})
    return @. aL[1]*ρN+aL[2]*τN-aL[3]*LapρN+aL[4]*ρN^(γ+1)+aL[5]*(ρN^2+2*ρp*ρn)
end

# Guleria Ver.
#function Calc_VΛΛ(aL,γ,ρN::Vector{Float64},ddρN::Vector{Float64},LapρN::Vector{Float64},τN::Vector{Float64},dτΛ::Vector{Float64})
#    return @. aL[1]*ρN+aL[2]*(ρN*dτΛ+τN)+aL[3]*(-LapρN+2*ddρN)+aL[4]*ρN^(γ+1)
#    #return @. aL[1]*ρN+aL[2]*(ρN*dτΛ+τN)-aL[3]*LapρN+aL[4]*ρN^(γ+1)
#end

function Calc_VΛN(aL,γ,ρN::Vector{Float64},ρΛ::Vector{Float64},LapρΛ::Vector{Float64},τΛ::Vector{Float64},ρq::Vector{Float64})
    return @. aL[1]*ρΛ+aL[2]*τΛ-aL[3]*LapρΛ+(γ+1)*aL[4]*(ρN^γ)*ρΛ+2*aL[5]*ρΛ*(ρN+ρq)
end

# Guleria Ver.
#function Calc_VΛN(aL,γ,ρΛ,τΛ,dτN,LapρΛ,ddρΛ,ρN)
#    return @. aL[1]*ρΛ+aL[2]*(τΛ+dτN*ρΛ)+aL[3]*(-LapρΛ+2*ddρΛ)+(γ+1)*aL[4]*(ρN^γ)*ρΛ
#    #return @. aL[1]*ρΛ+aL[2]*(τΛ+dτN*ρΛ)-aL[3]*LapρΛ+(γ+1)*aL[4]*(ρN^γ)*ρΛ
#end

function Calc_VNq(aN,σ,W0,ρN,ρq,τN,τq,LapρN,Lapρq,divJN,divJq)
    ans=zeros(Float64,Nmesh)
    @. ans[:]+=2*aN[1]*ρN[:] + 2*aN[2]*ρq[:] + (σ+2)*aN[3]*ρN[:]^(σ+1)
    @. ans[:]+=aN[4]*(σ*ρN[:]^(σ-1)*(ρq[:]^2+(ρN[:]-ρq[:])^2)+2*(ρN[:]^σ)*ρq[:])
    @. ans[:]+=aN[5]*τN[:]+aN[6]*τq[:]-2*aN[7]*LapρN[:]
    @. ans[:]+=-2*aN[8]*Lapρq[:]-0.5*W0*(divJN[:]+divJq[:])
    return ans
end

function Calc_Vcoul(ρp::Vector{Float64},rmesh,Z)
    Vcoul=zeros(Float64,Nmesh)
    Vcoul+=MyLib.SolvePoissonEq(ρp,rmesh,Z)
    @. Vcoul[2:Nmesh]=Vcoul[2:Nmesh]/rmesh[2:Nmesh]
    Vcoul[1]=InterPolEvenFunc0(Vcoul[2],Vcoul[3],Vcoul[4])
    #@. Vcoul[:]+=-(3*ρp[:]/π)^(1/3) #VB has no exchange term
    #Vcoul*=e2MeVfm/2 #Chabanat
    Vcoul*=e2MeVfm #Reainhard

    return Vcoul

end

function Calc_Wq(aN,W0,dρN::Vector{Float64},dρq::Vector{Float64},JN::Vector{Float64},Jq::Vector{Float64})
    return @. 0.5*W0*(dρN[:]+dρq[:])+2*aN[9]*JN[:]+2*aN[10]*Jq[:]
end


##########################################################
# HF iteration

function Calc_Density(Allocc,AllStates)
    ρ3=zeros(Float64,(3,Nmesh))
    dρ3=zeros(Float64,(3,Nmesh))
    Lapρ3=zeros(Float64,(3,Nmesh))
    τ3=zeros(Float64,(3,Nmesh))
    J3=zeros(Float64,(3,Nmesh))
    divJ3=zeros(Float64,(3,Nmesh))
    rmesh=getrmesh()

    for b in 1:3
        ρ3[b,:]=Calc_ρ(Allocc[b],AllStates[b],rmesh)
        dρ3[b,:]=Calc_dρ(ρ3[b,:],rmesh)
        Lapρ3[b,:]=Calc_Lapρ(ρ3[b,:],rmesh)
        τ3[b,:]=Calc_τ(Allocc[b],AllStates[b],rmesh)
        J3[b,:]=Calc_J(Allocc[b],AllStates[b],rmesh)
        divJ3[b,:]=Calc_divJ(J3[b,:],rmesh)
    end

    return ρ3,dρ3,Lapρ3,τ3,J3,divJ3
end

function Calc_Coef(ρ3,τ3,J3,aN,aL,pN,pΛ,AN::AtomNum)
	N=AN.N
	Z=AN.Z
	Λ=AN.Λ
    dρ3=zeros(Float64,(3,Nmesh))
    Lapρ3=zeros(Float64,(3,Nmesh))
    divJ3=zeros(Float64,(3,Nmesh))
    rmesh=getrmesh()

    for b in 1:3
        dρ3[b,:]=Calc_dρ(ρ3[b,:],rmesh)
        Lapρ3[b,:]=Calc_Lapρ(ρ3[b,:],rmesh)
        divJ3[b,:]=Calc_divJ(J3[b,:],rmesh)
    end

    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]

    dτ3=zeros(Float64,(3,Nmesh))
    ddρ3=zeros(Float64,(3,Nmesh))
    for b in 1:3
        dτ3[b,:]=Calc_dτ(τ3[b,:],rmesh)
        ddρ3[b,:]=Calc_ddρ(ρ3[b,:],rmesh)
    end
    dτN=dτ3[1,:]+dτ3[2,:]
    ddρN=ddρ3[1,:]+ddρ3[2,:]

    h2m=zeros(Float64,(3,Nmesh))
    dh2m=zeros(Float64,(3,Nmesh))
    ddh2m=zeros(Float64,(3,Nmesh))
    V=zeros(Float64,(3,Nmesh))
    W=zeros(Float64,(3,Nmesh))

    # Guleria
    #VΛΛ=Calc_VΛΛ(aL,pΛ.γ,ρN,ddρN,LapρN,τN,dτ3[3,:])
    #VΛN=Calc_VΛN(aL,pΛ.γ,ρ3[3,:],τ3[3,:],dτN,Lapρ3[3,:],ddρ3[3,:],ρN)
    # Rayet
    VΛΛ=Calc_VΛΛ(aL, pΛ.γ, ρN,LapρN,τN,ρ3[1,:],ρ3[2,:])
    VΛp=Calc_VΛN(aL, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[1,:])
    VΛn=Calc_VΛN(aL, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[2,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,Z)

    for b in 1:3
        if b==1 #proton
            h2m[b,:]+=Calc_h2mN(AN,b,aN,aL,ρN,ρ3[b,:],ρ3[3,:])
            dh2m[b,:]+=MyLib.diff1st5pt(h,h2m[b,:],1)
            ddh2m[b,:]+=MyLib.diff2nd5pt(h,h2m[b,:],1)
            V[b,:]+=VΛp+VNp+Vcoul
            W[b,:]+=Calc_Wq(aN,pN.W0,dρN,dρ3[b,:],JN,J3[b,:])
        elseif b==2 #neutron
            h2m[b,:]+=Calc_h2mN(AN,b,aN,aL,ρN,ρ3[b,:],ρ3[3,:])
            dh2m[b,:]+=MyLib.diff1st5pt(h,h2m[b,:],1)
            ddh2m[b,:]+=MyLib.diff2nd5pt(h,h2m[b,:],1)
            V[b,:]+=VΛn+VNn
            W[b,:]+=Calc_Wq(aN,pN.W0,dρN,dρ3[b,:],JN,J3[b,:])
        elseif b==3 #Lambda
            h2m[b,:]+=Calc_h2mΛ(AN,aL,ρN)
            dh2m[b,:]+=MyLib.diff1st5pt(h,h2m[b,:],1)
            ddh2m[b,:]+=MyLib.diff2nd5pt(h,h2m[b,:],1)
            V[b,:]+=VΛΛ
            @. W[b,:]+=0
        end
    end

    return h2m,dh2m,ddh2m,V,W
end

function CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh;rtol=1e-5)
    check=true
    diffρ=zeros(Float64,3)
    for b in 1:3
        ρold=Calc_ρ(Oldocc[b],OldStates[b],rmesh)
        ρnew=Calc_ρ(Newocc[b],NewStates[b],rmesh)
        #diffρ[b]=MyLib.IntTrap(rmesh,@. (4*π*rmesh[:]^2*(ρold[:]-ρnew[:])))
        diffρ[b]=MyLib.IntTrap(rmesh,(@. 4*π*rmesh[:]^2*(ρold[:]-ρnew[:])^2))
        if abs(diffρ[b])>rtol
            check=false
        end
    end

    println("diffρ=($(diffρ[1]),$(diffρ[2]),$(diffρ[3]))")

    return check

end

function HF_iter(AN::AtomNum;MaxIter=15,NParamType="SLy4",LParamType="HPL1",α=0.5)
    OldStates=InitialCondition(AN)
    Oldocc=Calc_occ(AN,OldStates)
    rmesh=getrmesh()
    Oldρ3,Olddρ3,OldLapρ3,Oldτ3,OldJ3,OlddivJ3=Calc_Density(Oldocc,OldStates)

    aN=NuclParameters.getaN(NParamType)
    aL=LambdaParameters.getaL(LParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(LParamType)

    #for debug
    #ρptest=zeros(Float64,(MaxIter,Nmesh))

    #plot()
    #calc several params
    for i in 1:MaxIter
        #for debug
        #ρptest[i,:]=Calc_ρ(Oldocc[1],OldStates[1],rmesh)

        h2m,dh2m,ddh2m,V,W=Calc_Coef(Oldρ3,Oldτ3,OldJ3,aN,aL,pN,pΛ,AN)

        NewStates=CalcAllStates(h2m,dh2m,ddh2m,V,W,rmesh)
        Newocc=Calc_occ(AN,NewStates)
        Newρ3,Newdρ3,NewLapρ3,Newτ3,NewJ3,NewdivJ3=Calc_Density(Newocc,NewStates)

        if CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh)==true
            return Newocc,NewStates
            break
        end
        println(i)
        OldStates=NewStates
        Oldocc=Newocc
        Oldρ3=Oldρ3*(1-α)+Newρ3*α
        Oldτ3=Oldτ3*(1-α)+Newτ3*α
        OldJ3=OldJ3*(1-α)+NewJ3*α
    end

    #=
    plot(xlabel="r",ylabel="ρn")
    for i in 1:MaxIter
        plot!(rmesh,ρptest[i,:],label="$i")
    end
    plot!()
    =#

end

############################################3
# out put files

function OutPutFiles(AN::AtomNum;NParamType="SLy4",LParamType="HPL1", α=0.5)
    Ansocc,AnsStates=HF_iter(AN,NParamType=NParamType,LParamType=LParamType,MaxIter=50,α=α)

    Z=AN.Z
    N=AN.N
    Λ=AN.Λ
    rm("data/Z$(Z)N$(N)L$(Λ)_$(NParamType)$(LParamType)",force=true,recursive=true)
    mkpath("data/Z$(Z)N$(N)L$(Λ)_$(NParamType)$(LParamType)")
    cd("data/Z$(Z)N$(N)L$(Λ)_$(NParamType)$(LParamType)")
    WriteStates(AN,Ansocc,AnsStates,NParamType,LParamType)
    WriteWaveFunc(AN,Ansocc,AnsStates,NParamType,LParamType)
    WriteDensityPot(AN,Ansocc,AnsStates,NParamType,LParamType)
	if Λ==1
        WriteTotalEnergy(AN,Ansocc,AnsStates,NParamType,LParamType)
	elseif Λ==0
		WriteTotalEnergy(AN,Ansocc,AnsStates,NParamType)
	end
    cd("../..")
end

function WriteHeader(io::IOStream,AN,NParamType,LParamType)
	rmesh=getrmesh()
    Z=AN.Z
    N=AN.N
    Λ=AN.Λ
    write(io, "# Nuclear Parameter=$(NParamType)\n")
    write(io, "# Lambda Parameter=$(LParamType)\n")
    write(io, "# Z=$(Z), N=$(N), Λ=$(Λ)\n")
    write(io, "# Number of mesh=$(Nmesh)\n")
    write(io, "# rmax=$(rmax)\n")
    write(io, "# Matching point of shooting = $(rmesh[Nmatch])\n\n")
end

function WriteHeader(io::IOStream,AN,NParamType)
	rmesh=getrmesh()
    Z=AN.Z
    N=AN.N
	Λ=AN.Λ
    write(io, "# Nuclear Parameter=$(NParamType)\n")
    write(io, "# Z=$(Z), N=$(N), Λ=$(Λ)\n")
    write(io, "# Number of mesh=$(Nmesh)\n")
    write(io, "# rmax=$(rmax)\n")
    write(io, "# Matching point of shooting = $(rmesh[Nmatch])\n\n")
end

function WriteStates(AN::AtomNum,Ansocc,AnsStates,NParamType,LParamType)
    io=open("states.csv","w")
	WriteHeader(io,AN,NParamType,LParamType)

	write(io, "Baryon Type, occ, j, l, Energy(MeV)\n")
    for b=1:3
        for i=eachindex(AnsStates[b])
            j=AnsStates[b][i].QN.j
            l=AnsStates[b][i].QN.l
            E=AnsStates[b][i].E
            if b==1 #proton
                write(io, "proton, $(Ansocc[b][i]), $(j), $(l), $(E)\n")
            elseif b==2 #neutron
                write(io, "neutron, $(Ansocc[b][i]), $(j), $(l), $(E)\n")
            elseif b==3 #Λ
                write(io, "lambda, $(Ansocc[b][i]), $(j), $(l), $(E)\n")
            end
        end
    end
    close(io)
end

function WriteWaveFunc(AN,Ansocc,AnsStates,NParamType,LParamType)
    io=open("wavefunc.csv","w")
	WriteHeader(io,AN,NParamType,LParamType)
	rmesh=getrmesh()

    write(io, "r(fm)")
    for b in 1:3
        for i=eachindex(AnsStates[b])
            if b==1
                write(io, ",Rp$(i)")
            elseif b==2
                write(io, ",Rn$(i)")
            elseif b==3
                write(io, ",Rl$(i)")
            end
        end
    end

    write(io, "\n")
    for n in 1:Nmesh
        write(io, "$(rmesh[n])")
        for b in 1:3
            for i=eachindex(AnsStates[b])
                write(io, ",$(AnsStates[b][i].ψ[n])")
            end
        end
        write(io, "\n")
    end

    close(io)

end

function WriteDensityPot(AN,Ansocc,AnsStates,NParamType,LParamType)
    io1=open("density.csv","w")
    WriteHeader(io1,AN,NParamType,LParamType)
	rmesh=getrmesh()

    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Ansocc,AnsStates)
    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]

    aN=NuclParameters.getaN(NParamType)
    aL=LambdaParameters.getaL(LParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(LParamType)
    VΛΛ=Calc_VΛΛ(aL, pΛ.γ, ρN,LapρN,τN,ρ3[1,:],ρ3[2,:])
    VΛp=Calc_VΛN(aL, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[1,:])
    VΛn=Calc_VΛN(aL, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[2,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,AN.Z)

    write(io1, "r(fm)")
    write(io1, ",Rhop,dRhop,LapRhop,Jp,DivJp")
    write(io1, ",Rhon,dRhon,LapRhon,Jn,DivJn")
    write(io1, ",Rhol,dRhol,LapRhol,Jl,DivJl")
    write(io1, ",RhoN,dRhoN,LapRhoN,JN,DivJN\n")

    for n in 1:Nmesh
        write(io1, "$(rmesh[n])")
        for b in 1:4
            if b<=3
                write(io1, ",$(ρ3[b,n])")
                write(io1, ",$(dρ3[b,n])")
                write(io1, ",$(Lapρ3[b,n])")
                write(io1, ",$(J3[b,n])")
                write(io1, ",$(divJ3[b,n])")
            elseif b==4
                write(io1, ",$(ρN[n])")
                write(io1, ",$(dρN[n])")
                write(io1, ",$(LapρN[n])")
                write(io1, ",$(JN[n])")
                write(io1, ",$(divJN[n])")
            end
        end
        write(io1, "\n")
    end

    close(io1)

    io2=open("potential.csv","w")
    WriteHeader(io2,AN,NParamType,LParamType)

    write(io2, "r(fm)")
    write(io2, ",Vll(MeV),Vlp(MeV),Vln(MeV),VNp(MeV),VNn(MeV),Vcoul(MeV)\n")
    for n in 1:Nmesh
        write(io2, "$(rmesh[n])")
        write(io2, ",$(VΛΛ[n])")
        write(io2, ",$(VΛp[n])")
        write(io2, ",$(VΛn[n])")
        write(io2, ",$(VNp[n])")
        write(io2, ",$(VNn[n])")
        write(io2, ",$(Vcoul[n])")
        write(io2, "\n")
    end

    close(io2)

end

##################################################
# Calculate Binding Energy
function Hamiltonian_N(aN,σ,W0,ρ3,ρN,τ3,τN,Lapρ3,LapρN,J3,JN,divJ3,divJN)
    Hn=zeros(Float64,Nmesh)
    @. Hn += ħc^2/(2*mpMeV)*τ3[1,:] + ħc^2/(2*mnMeV)*τ3[2,:]
    @. Hn += aN[1]*ρN[:]^2
    @. Hn += aN[2]*(ρ3[1,:]^2 + ρ3[2,:]^2)
    @. Hn += aN[3]*ρN[:]^(σ+2)
    @. Hn += aN[4]*ρN[:]^σ*(ρ3[1,:]^2 + ρ3[2,:]^2)
    @. Hn += aN[5]*τN[:]*ρN[:]
    @. Hn += aN[6]*(τ3[1,:]*ρ3[1,:] + τ3[2,:]*ρ3[2,:])
    @. Hn += -aN[7]*ρN[:]*LapρN[:]
    @. Hn += -aN[8]*(ρ3[1,:]*Lapρ3[1,:] + ρ3[2,:]*Lapρ3[2,:])
	@. Hn += -0.5*W0*(divJN[:]*ρN[:] + divJ3[1,:]*ρ3[1,:] + divJ3[2,:]*ρ3[2,:])
    @. Hn += aN[9]*JN[:]^2
    @. Hn += aN[10]*(J3[1,:]^2 + J3[2,:]^2)
    return Hn
end

function Energy_N(aN,σ,W0,ρ3,ρN,τ3,τN,Lapρ3,LapρN,J3,JN,divJ3,divJN)
	rmesh=getrmesh()
    Hn=Hamiltonian_N(aN,σ,W0,ρ3,ρN,τ3,τN,Lapρ3,LapρN,J3,JN,divJ3,divJN)
    En=MyLib.IntTrap(rmesh,(@. rmesh[:]^2*Hn[:]))*4*π
    return En
end

function Hamiltonian_L(aL,γ,ρ3,ρN,τ3,τN,Lapρ3,LapρN)
    Hl=zeros(Float64,Nmesh)
    @. Hl += ħc^2/(2*mΛMeV)*τ3[3,:]
    @. Hl += aL[1]*ρN[:]*ρ3[3,:]
    @. Hl += aL[2]*(τ3[3,:]*ρN[:] + τN*ρ3[3,:])
    @. Hl -= aL[3]*(ρ3[3,:]*LapρN[:])
    @. Hl += aL[4]*ρN[:]^(γ+1)*ρ3[3,:]
    @. Hl += aL[5]*ρ3[3,:]*(ρN[:]^2 + 2*ρ3[1,:]*ρ3[2,:])
    return Hl
end

function Energy_L(aL,γ,ρ3,ρN,τ3,τN,Lapρ3,LapρN)
	rmesh=getrmesh()
    Hl=Hamiltonian_L(aL,γ,ρ3,ρN,τ3,τN,Lapρ3,LapρN)
    En=MyLib.IntTrap(rmesh,(@. rmesh[:]^2*Hl[:]))*4*π
    return En
end

function H_coul_dir(ρp,rmesh,Z)
    Hcoul_dir=zeros(Float64,Nmesh)
    Hcoul_dir+=MyLib.SolvePoissonEq(ρp,rmesh,Z)
    @. Hcoul_dir[2:Nmesh]*=0.5*ρp[2:Nmesh]/rmesh[2:Nmesh]
    Hcoul_dir[1]=InterPolEvenFunc0(Hcoul_dir[2],Hcoul_dir[3],Hcoul_dir[4])
    #@. Hcoul_exch[:] -= 0.75*ρp[:]*(3*ρp[:]/π)^(1/3)
    #Hcoul_dir*=e2MeVfm/2 #Chabanatd
    Hcoul_dir*=e2MeVfm #Reainhard

    return Hcoul_dir
end

function Energy_coul_dir(ρp,rmesh,Z)
	Hc_dir=H_coul_dir(ρp,rmesh,Z)
	Ec_dir=MyLib.IntTrap(rmesh,(@. rmesh[:]^2*Hc_dir[:]))*4*π
	return Ec_dir
end

function H_coul_exch(ρp)
    Hcoul_exch=zeros(Float64,Nmesh)
    @. Hcoul_exch[:] -= 0.75*ρp[:]*(3*ρp[:]/π)^(1/3)
    #Hcoul_exch*=e2MeVfm/2 #Chabanatd
    Hcoul_exch*=e2MeVfm #Reainhard

    return Hcoul_exch
end

function Energy_coul_exch(ρp)
	rmesh=getrmesh()
	Hc_exch=H_coul_exch(ρp)
	Ec_dir=MyLib.IntTrap(rmesh,(@. rmesh[:]^2*Hc_exch[:]))*4*π
	return Ec_dir
end

function Energy_Pair()
    Ep=0.0
    return Ep
end

function Energy_CM_dir(τ3)
	rmesh=getrmesh()
    Ecm_dir=0.0
    Z=AN.Z
    N=AN.N
    Λ=AN.Λ
    for b in 1:3
        Ecm_dir += MyLib.IntTrap(rmesh,@. rmesh[:]^2*τ3[b,:])*4*π
    end
    Ecm_dir*=ħc^2/(2*(mpMeV*Z + mnMeV*N + mΛMeV*Λ))
    return Ecm_dir
end

function Energy_CM_exch()
    Ecm_exch=0.0
    return Ecm_exch
end

function Energy_Kin(AN,τ3)
	E_Kin=0.0
	rmesh=getrmesh()
	Z=AN.Z
	N=AN.N
	for b in 1:3
		QN=QuantumNumber(0.5,0,b)
		mass=getmass(QN)
		E_Kin += MyLib.IntTrap(rmesh,@. rmesh[:]^2*(ħc^2/(2*mass)-ħc^2/(2*(mpMeV*Z+mnMeV*N)))*τ3[b,:] )*4*π
	end
	return E_Kin
end

function Energy_N_R(aN,σ,ρ3,ρN)
	Hn_R=zeros(Float64,Nmesh)
	rmesh=getrmesh()
	@. Hn_R += aN[3]*ρN[:]^(σ+2)
    @. Hn_R += aN[4]*ρN[:]^σ*(ρ3[1,:]^2 + ρ3[2,:]^2)
	En_R=0.5*σ*MyLib.IntTrap(rmesh,@. rmesh[:]^2*Hn_R[:])*4*π

	return En_R
end

function Energy_L_R()

end

function Energy_SPS(Ansocc,AnsStates)
	E=0.0
	for b in 1:3
		for i=eachindex(Ansocc[b])
			j=AnsStates[b][i].QN.j
			E+=Ansocc[b][i]*AnsStates[b][i].E*(2*j+1)
		end
	end
	return E
end

function WriteTotalEnergy(AN,Ansocc,AnsStates,NParamType)
    io1=open("Energy1.csv","w")
    WriteHeader(io1,AN,NParamType)
	rmesh=getrmesh()

    aN=NuclParameters.getaN(NParamType)
    pN=NuclParameters.getParams(NParamType)

    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Ansocc,AnsStates)
    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]

	write(io1,"Etot/A(MeV), Etot1(MeV), 	EN(MeV),	 Ec_dir(MeV),	Ecm_dir(MeV)\n")

	#directory integrate energy density functional
	En=Energy_N(aN,pN.σ,pN.W0,ρ3,ρN,τ3,τN,Lapρ3,LapρN,J3,JN,divJ3,divJN)
	Ec_dir=Energy_coul_dir(ρ3[1,:],rmesh,AN.Z)
	#Ec_exch=Energy_coul_exch(ρ3[1,:])
	#Epair=Energy_Pair()
	Ecm_dir=Energy_CM_dir(τ3)
	#Ecm_exch=Energy_CM_exch()
	Etot1=En + Ec_dir - Ecm_dir #+ Ec_exch + Epair  - Ecm_exch

	write(io1,"$(Etot1/(AN.N+AN.Z))")
	write(io1,",$(Etot1)")
	write(io1,",$(En)")
	write(io1,",$(Ec_dir)")
	#write(io1,",$(Ec_exch)")
	#write(io1,",$(Epair)")
	write(io1,",$(Ecm_dir)")
	#write(io1,",$(Ecm_exch)\n")
	close(io1)

	#using single-particle energy and rearrangement energy
	io2=open("Energy2.csv","w")
    WriteHeader(io2,AN,NParamType)

	E_Kin=Energy_Kin(AN,τ3)
	E_SPS=Energy_SPS(Ansocc,AnsStates)
	En_R=Energy_N_R(aN,pN.σ,ρ3,ρN)
	Etot2=0.5*(E_Kin+E_SPS)- En_R# + Epair - Ecm_dir - Ecm_exch

	write(io2, "Etot2/A(MeV), Etot2(MeV), EKin(MeV), ESPS(MeV), ER(MeV)\n")
	write(io2,"$(Etot2/(AN.N+AN.Z))")
	write(io2,",$(Etot2)")
	write(io2,",$(E_Kin)")
	write(io2,",$(E_SPS)")
	write(io2,",$(En_R)\n")

	close(io2)

end