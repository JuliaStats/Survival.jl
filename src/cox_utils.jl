# Compute first and last
firsts(s) = [s[t].status && (t==1 || s[t] > s[t-1]) for t = 1:length(s)]
lasts(s) = [s[t].status && (t==length(s) || s[t+1] > s[t]) for t = 1:length(s)]


# Compute "complementary cumulative sum" on first dimension
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
