module LinearElasticity

using ConstructionBase: constructorof
using Tensorial: SymmetricSecondOrderTensor, SymmetricFourthOrderTensor, Vec

export TensorStress,
    TensorStrain,
    StiffnessTensor,
    ComplianceTensor,
    EngineeringStress,
    EngineeringStrain,
    ComplianceMatrix,
    StiffnessMatrix,
    Cubic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Orthorhombic,
    Monoclinic,
    Triclinic

"Represent one of the seven crystal systems."
abstract type CrystalSystem end
"""
    Triclinic()
Represent the triclinic system.
"""
struct Triclinic <: CrystalSystem end
"""
    Monoclinic()
Represent the monoclinic system.
"""
struct Monoclinic <: CrystalSystem end
"""
    Orthorhombic()
Represent the orthorhombic system.
"""
struct Orthorhombic <: CrystalSystem end
"""
    Tetragonal()
Represent the tetragonal system.
"""
struct Tetragonal <: CrystalSystem end
"""
    Cubic()
Represent the cubic system.
"""
struct Cubic <: CrystalSystem end
"""
    Trigonal()
Represent the trigonal system.
"""
struct Trigonal <: CrystalSystem end
"""
    Hexagonal()
Represent the hexagonal system.
"""
struct Hexagonal <: CrystalSystem end

abstract type Stress{T,N} <: AbstractArray{T,N} end
abstract type Strain{T,N} <: AbstractArray{T,N} end
abstract type Stiffness{T,N} <: AbstractArray{T,N} end
abstract type Compliance{T,N} <: AbstractArray{T,N} end
struct TensorStress{T} <: Stress{T,2}
    data::SymmetricSecondOrderTensor{3,T,6}
end
TensorStress(m::AbstractMatrix) = TensorStress(SymmetricSecondOrderTensor{3}(m))
TensorStress(data...) = TensorStress(SymmetricSecondOrderTensor{3}(data...))
struct TensorStrain{T} <: Strain{T,2}
    data::SymmetricSecondOrderTensor{3,T,6}
end
TensorStrain(m::AbstractMatrix) = TensorStrain(SymmetricSecondOrderTensor{3}(m))
TensorStrain(data...) = TensorStrain(SymmetricSecondOrderTensor{3}(data...))
struct StiffnessTensor{T} <: Stiffness{T,4}
    data::SymmetricFourthOrderTensor{3,T}
end
struct ComplianceTensor{T} <: Compliance{T,4}
    data::SymmetricFourthOrderTensor{3,T}
end
struct EngineeringStress{T} <: Stress{T,1}
    data::Vec{6,T}
end
EngineeringStress(x::TensorStress) = convert(EngineeringStress{eltype(x)}, x)
EngineeringStress(v::AbstractVector) = EngineeringStress(Vec{6}(v))
EngineeringStress(data...) = EngineeringStress(Vec{6}(data...))
TensorStress(x::EngineeringStress) = convert(TensorStress{eltype(x)}, x)
struct EngineeringStrain{T} <: Strain{T,1}
    data::Vec{6,T}
end
EngineeringStrain(x::TensorStrain) = convert(EngineeringStrain{eltype(x)}, x)
EngineeringStrain(v::AbstractVector) = EngineeringStrain(Vec{6}(v))
EngineeringStrain(data...) = EngineeringStrain(Vec{6}(data...))
TensorStrain(x::EngineeringStrain) = convert(TensorStrain{eltype(x)}, x)
struct StiffnessMatrix{T} <: Stiffness{T,2}
    data::SymmetricSecondOrderTensor{6,T,21}
end
StiffnessMatrix(x::StiffnessTensor) = convert(StiffnessMatrix{eltype(x)}, x)
StiffnessMatrix(m::AbstractMatrix) = StiffnessMatrix(SymmetricSecondOrderTensor{6}(m))
StiffnessMatrix(data...) = StiffnessMatrix(SymmetricSecondOrderTensor{6}(data...))
StiffnessTensor(x::StiffnessMatrix) = convert(StiffnessTensor{eltype(x)}, x)
struct ComplianceMatrix{T} <: Compliance{T,2}
    data::SymmetricSecondOrderTensor{6,T,21}
end
ComplianceMatrix(x::ComplianceTensor) = convert(ComplianceMatrix{eltype(x)}, x)
ComplianceMatrix(m::AbstractMatrix) = ComplianceMatrix(SymmetricSecondOrderTensor{6}(m))
ComplianceMatrix(data...) = ComplianceMatrix(SymmetricSecondOrderTensor{6}(data...))
ComplianceTensor(x::ComplianceMatrix) = convert(ComplianceTensor{eltype(x)}, x)

Base.size(::Union{TensorStress,TensorStrain}) = (3, 3)
Base.size(::Union{StiffnessTensor,ComplianceTensor}) = (3, 3, 3, 3)
Base.size(::Union{EngineeringStress,EngineeringStrain}) = (6,)
Base.size(::Union{StiffnessMatrix,ComplianceMatrix}) = (6, 6)

Base.getindex(A::Union{Stress,Strain,Stiffness,Compliance}, i...) = getindex(A.data, i...)

Base.parent(A::Union{Stress,Strain,Stiffness,Compliance}) = A.data

Base.IndexStyle(::Type{<:Union{Stress,Strain,Stiffness,Compliance}}) = IndexLinear()

# See https://github.com/JuliaLang/julia/blob/cb9acf5/base/arraymath.jl#L19-L26
for op in (:*, :/)
    @eval Base.$op(A::Union{Stress,Strain,Stiffness,Compliance}, B::Number) =
        constructorof(typeof(A))(Base.broadcast_preserving_zero_d($op, A, B))
end
Base.:*(B::Number, A::Union{Stress,Strain,Stiffness,Compliance}) = A * B

Base.:-(A::Union{Stress,Strain,Stiffness,Compliance}) =
    constructorof(typeof(A))(Base.broadcast_preserving_zero_d(-, A))

for op in (:+, :-)
    for T in (
        :TensorStress,
        :TensorStrain,
        :EngineeringStress,
        :EngineeringStrain,
        :StiffnessMatrix,
        :ComplianceMatrix,
        :StiffnessTensor,
        :ComplianceTensor,
    )
        @eval Base.$op(A::$T, B::$T) = $T(Base.broadcast_preserving_zero_d($op, A, B))
    end
end

include("conversion.jl")
include("invariants.jl")
include("SymmetryCriteria.jl")
# include("StabilityCriteria.jl")
include("Isotropic.jl")
include("misc.jl")
include("energy.jl")
include("solve.jl")
include("ULICS.jl")

end
