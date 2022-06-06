using Compat.LinearAlgebra: dot

export elastic_energy_density

elastic_energy_density(σ::EngineeringStress, ϵ::EngineeringStrain) = dot(σ, ϵ) / 2
elastic_energy_density(σ::TensorStress, ε::TensorStrain) = double_contraction(σ, ε) / 2
elastic_energy_density(ε, σ) = elastic_energy_density(ε, σ)
elastic_energy_density(cᵢⱼ::StiffnessMatrix, ϵ::EngineeringStrain) = dot(ϵ, cᵢⱼ, ϵ) / 2
elastic_energy_density(sᵢⱼ::ComplianceMatrix, σ::EngineeringStress) = dot(σ, sᵢⱼ, σ) / 2
elastic_energy_density(cᵢⱼₖₗ::StiffnessTensor, ε::TensorStrain) =
    elastic_energy_density(StiffnessMatrix(cᵢⱼₖₗ), EngineeringStrain(ε))
elastic_energy_density(sᵢⱼₖₗ::ComplianceTensor, σ::TensorStress) =
    elastic_energy_density(ComplianceMatrix(sᵢⱼₖₗ), EngineeringStress(σ))
