

## Tables.jl integration??
## get it basically for free....
Tables.istable(::Type{<:NonparametricEstimator}) = true
Tables.columnaccess(::Type{<:NonparametricEstimator}) = true
Tables.columns(np::T) where {T<:NonparametricEstimator} = (;Iterators.zip(propertynames(np),map(nm->getproperty(np,nm),propertynames(np)))...)
Tables.schema(np::T) where {T<:NonparametricEstimator} = Tables.Schema(propertynames(np),eltype.(fieldtypes(typeof(np))))


# km = fit(KaplanMeier,ev) |> DataFrame