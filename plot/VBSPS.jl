mutable struct SPS
	τ::Int64 #τ=1::neutron, τ=-1::protons
	l::Int64
	j::Float64
	E::Float64
end

O_SPS=[
	SPS(-1,0,0.5,28.76)
]