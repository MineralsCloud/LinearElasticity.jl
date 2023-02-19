module LinearElasticity

using LinearElasticityBase:
    Stress,
    Strain,
    Stiffness,
    Compliance,
    TensorStress,
    TensorStrain,
    StiffnessTensor,
    ComplianceTensor,
    EngineeringStress,
    EngineeringStrain,
    ComplianceMatrix,
    StiffnessMatrix,
    energydensity,
    principal_values,
    principal_axes,
    principal_invariants,
    main_invariants,
    hydrostatic,
    deviatoric,
    contraction,
    double_contraction,
    rotate_axes,
    rotate,
    to_tensor,
    to_voigt,
    isequivalent,
    isbiaxial,
    isuniaxial,
    ⊡,
    ⩵

export Stress,
    Stress,
    Strain,
    Stiffness,
    Compliance,
    TensorStress,
    TensorStrain,
    StiffnessTensor,
    ComplianceTensor,
    EngineeringStress,
    EngineeringStrain,
    ComplianceMatrix,
    StiffnessMatrix,
    energydensity,
    principal_values,
    principal_axes,
    principal_invariants,
    main_invariants,
    hydrostatic,
    deviatoric,
    contraction,
    double_contraction,
    rotate_axes,
    rotate,
    to_tensor,
    to_voigt,
    isequivalent,
    isbiaxial,
    isuniaxial,
    ⊡,
    ⩵
export CubicConstraint,
    HexagonalConstraint,
    TetragonalConstraint,
    TrigonalConstraint,
    OrthorhombicConstraint,
    MonoclinicConstraint,
    TriclinicConstraint

abstract type SymmetryConstraint end
struct TriclinicConstraint <: SymmetryConstraint end
struct MonoclinicConstraint <: SymmetryConstraint end
struct OrthorhombicConstraint <: SymmetryConstraint end
struct TetragonalConstraint <: SymmetryConstraint end
struct CubicConstraint <: SymmetryConstraint end
struct TrigonalConstraint <: SymmetryConstraint end
struct HexagonalConstraint <: SymmetryConstraint end

include("SymmetryCriteria.jl")
# include("StabilityCriteria.jl")
include("Isotropic.jl")
include("misc.jl")
include("Solve.jl")
include("ULICS.jl")

end
