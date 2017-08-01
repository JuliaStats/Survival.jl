abstract type SurvivalModel end

function Base.show(io::IO, obj::SurvivalModel)
    print(io,"\nModel: ", obj.model, obj.formula,"\n\n")
    print(io,obj.coefmat)
end

StatsBase.coef(SM::SurvivalModel) = SM.coefmat.cols[1]

struct CoxModel <: SurvivalModel
    model::AbstractString
    formula::Formula
    coefmat::CoefTable
    M::ModelFrame
    loglik::Float64
    score::Array{Float64,1}
    fischer_info::Array{Float64,2}
end
