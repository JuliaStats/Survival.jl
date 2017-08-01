# Compute first and last! Could be optimized!
firsts(S) = [S[t].status && (t==1 || S[t] > S[t-1]) for t = 1:length(S)]
lasts(S) = [S[t].status && (t==length(S) || S[t+1] > S[t]) for t = 1:length(S)]


# Computes "complementary cumulative sum" on first dimension
function after(v)
    afterv = init_after(v)
    after!(afterv,v)
    return afterv
end


function init_after(v)
    newsize = collect(size(v))
    newsize[1] += 1
    afterv = zeros(eltype(v),newsize...)
    return afterv
end

function after!{N,T}(afterv, v::AbstractArray{T,N})
    trails = fill(:,N-1)
    cumsum!(@view(afterv[(end-1):-1:1,trails...]),@view(v[end:-1:1,trails...]),1)
    return
end
