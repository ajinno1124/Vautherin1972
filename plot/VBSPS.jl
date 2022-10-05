mutable struct SPS
	τ::String
	l::Int64
	j::Float64
	E::Float64
end

O_SPS=[
	SPS(proton,0,0.5,28.76),
	SPS(proton,1,1.5,16.77),
	SPS(proton,1,0.5,12.97),
	SPS(proton,2,2.5,1.78),
	SPS(neutron,0,0.5,32.96),
	SPS(neutron,1,1.5,20.81),
	SPS()
]