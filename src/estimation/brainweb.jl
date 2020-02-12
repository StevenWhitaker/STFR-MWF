"""
    BrainWeb

Abstract type for storing BrainWeb tissue parameters.
"""
abstract type BrainWeb{N} end

"""
    BrainWeb3Comp <: BrainWeb

Type in which to store three-compartment tissue data using the BrainWeb phantom.
"""
struct BrainWeb3Comp{N} <: BrainWeb{N}
    tissue::Array{Symbol,N}
    M0::Array{Float64,N}
    ff::Array{Float64,N}
    fm::Array{Float64,N}
    T1f::Array{Float64,N}
    T1s::Array{Float64,N}
    T1m::Array{Float64,N}
    T2f::Array{Float64,N}
    T2s::Array{Float64,N}
    T2m::Array{Float64,N}
    τfs::Array{Float64,N}
    τfm::Array{Float64,N}
    Δωf::Array{Float64,N}

    BrainWeb3Comp(id::Array{UInt8,N}) where {N} = begin
        new{N}(
            Array{Symbol,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id))
        )
    end
end

"""
    BrainWeb4Comp <: BrainWeb

Type in which to store four-compartment tissue data using the BrainWeb phantom.
"""
struct BrainWeb4Comp{N} <: BrainWeb{N}
    tissue::Array{Symbol,N}
    M0::Array{Float64,N}
    ff::Array{Float64,N}
    fm::Array{Float64,N}
    fa::Array{Float64,N}
    T1f::Array{Float64,N}
    T1s::Array{Float64,N}
    T1m::Array{Float64,N}
    T1a::Array{Float64,N}
    T2f::Array{Float64,N}
    T2s::Array{Float64,N}
    T2m::Array{Float64,N}
    T2a::Array{Float64,N}
    τfs::Array{Float64,N}
    τfm::Array{Float64,N}
    τfa::Array{Float64,N}
    Δωf::Array{Float64,N}

    BrainWeb4Comp(id::Array{UInt8,N}) where {N} = begin
        new{N}(
            Array{Symbol,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id))
        )
    end
end

"""
    BrainWeb9Comp <: BrainWeb

Type in which to store nine-compartment tissue data using the BrainWeb phantom.
"""
struct BrainWeb9Comp{N} <: BrainWeb{N}
    tissue::Array{Symbol,N}
    M0::Array{Float64,N}
    ff1::Array{Float64,N}
    ff2::Array{Float64,N}
    ff3::Array{Float64,N}
    fs1::Array{Float64,N}
    fs2::Array{Float64,N}
    fs3::Array{Float64,N}
    fm1::Array{Float64,N}
    fm2::Array{Float64,N}
    fm3::Array{Float64,N}
    T1f1::Array{Float64,N}
    T1f2::Array{Float64,N}
    T1f3::Array{Float64,N}
    T1s1::Array{Float64,N}
    T1s2::Array{Float64,N}
    T1s3::Array{Float64,N}
    T1m1::Array{Float64,N}
    T1m2::Array{Float64,N}
    T1m3::Array{Float64,N}
    T2f1::Array{Float64,N}
    T2f2::Array{Float64,N}
    T2f3::Array{Float64,N}
    T2s1::Array{Float64,N}
    T2s2::Array{Float64,N}
    T2s3::Array{Float64,N}
    T2m1::Array{Float64,N}
    T2m2::Array{Float64,N}
    T2m3::Array{Float64,N}
    τf1f2::Array{Float64,N}
    τf1f3::Array{Float64,N}
    τf1s1::Array{Float64,N}
    τf1s2::Array{Float64,N}
    τf1s3::Array{Float64,N}
    τf1m1::Array{Float64,N}
    τf1m2::Array{Float64,N}
    τf1m3::Array{Float64,N}
    τf2f3::Array{Float64,N}
    τf2s1::Array{Float64,N}
    τf2s2::Array{Float64,N}
    τf2s3::Array{Float64,N}
    τf2m1::Array{Float64,N}
    τf2m2::Array{Float64,N}
    τf2m3::Array{Float64,N}
    τf3s1::Array{Float64,N}
    τf3s2::Array{Float64,N}
    τf3s3::Array{Float64,N}
    τf3m1::Array{Float64,N}
    τf3m2::Array{Float64,N}
    τf3m3::Array{Float64,N}
    τs1s2::Array{Float64,N}
    τs1s3::Array{Float64,N}
    τs2s3::Array{Float64,N}
    τm1m2::Array{Float64,N}
    τm1m3::Array{Float64,N}
    τm2m3::Array{Float64,N}
    Δωf::Array{Float64,N}

    BrainWeb9Comp(id::Array{UInt8,N}) where {N} = begin
        new{N}(
            Array{Symbol,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id)),
            Array{Float64,N}(undef, size(id))
        )
    end
end

Base.getproperty(obj::BrainWeb9Comp, s::Symbol) = begin
    if s == :ff
        return sum(map(s -> getfield(obj, s), [:ff1, :ff2, :ff3]))
    else
        return getfield(obj, s)
    end
end

Base.size(obj::BrainWeb) = size(obj.tissue)

function getmask(obj::BrainWeb{N}, tissues::Symbol...) where {N}

    mask = dropdims(any(cat([obj.tissue .== t for t in tissues]..., dims = N+1), dims = N+1), dims = N+1)

end

function getparams(obj::BrainWeb3Comp{N}, tissues::Symbol...) where {N}

    mask = getmask(obj, tissues...)
    return @views((mask,
                   obj.M0[mask],
                   obj.ff[mask],
                   obj.fm[mask],
                   obj.T1f[mask],
                   obj.T1s[mask],
                   obj.T1m[mask],
                   obj.T2f[mask],
                   obj.T2s[mask],
                   obj.T2m[mask],
                   obj.τfs[mask],
                   obj.τfm[mask],
                   obj.Δωf[mask]))

end

function getparams(obj::BrainWeb4Comp{N}, tissues::Symbol...) where {N}

    mask = getmask(obj, tissues...)
    return @views((mask,
                   obj.M0[mask],
                   obj.ff[mask],
                   obj.fm[mask],
                   obj.fa[mask],
                   obj.T1f[mask],
                   obj.T1s[mask],
                   obj.T1m[mask],
                   obj.T1a[mask],
                   obj.T2f[mask],
                   obj.T2s[mask],
                   obj.T2m[mask],
                   obj.T2a[mask],
                   obj.τfs[mask],
                   obj.τfm[mask],
                   obj.τfa[mask],
                   obj.Δωf[mask]))

end

function getparams(obj::BrainWeb9Comp{N}, tissues::Symbol...) where {N}

    mask = getmask(obj, tissues...)
    return @views((mask,
                   obj.M0[mask],
                   obj.ff1[mask],
                   obj.ff2[mask],
                   obj.ff3[mask],
                   obj.fs1[mask],
                   obj.fs2[mask],
                   obj.fs3[mask],
                   obj.fm1[mask],
                   obj.fm2[mask],
                   obj.fm3[mask],
                   obj.T1f1[mask],
                   obj.T1f2[mask],
                   obj.T1f3[mask],
                   obj.T1s1[mask],
                   obj.T1s2[mask],
                   obj.T1s3[mask],
                   obj.T1m1[mask],
                   obj.T1m2[mask],
                   obj.T1m3[mask],
                   obj.T2f1[mask],
                   obj.T2f2[mask],
                   obj.T2f3[mask],
                   obj.T2s1[mask],
                   obj.T2s2[mask],
                   obj.T2s3[mask],
                   obj.T2m1[mask],
                   obj.T2m2[mask],
                   obj.T2m3[mask],
                   obj.τf1f2[mask],
                   obj.τf1f3[mask],
                   obj.τf1s1[mask],
                   obj.τf1s2[mask],
                   obj.τf1s3[mask],
                   obj.τf1m1[mask],
                   obj.τf1m2[mask],
                   obj.τf1m3[mask],
                   obj.τf2f3[mask],
                   obj.τf2s1[mask],
                   obj.τf2s2[mask],
                   obj.τf2s3[mask],
                   obj.τf2m1[mask],
                   obj.τf2m2[mask],
                   obj.τf2m3[mask],
                   obj.τf3s1[mask],
                   obj.τf3s2[mask],
                   obj.τf3s3[mask],
                   obj.τf3m1[mask],
                   obj.τf3m2[mask],
                   obj.τf3m3[mask],
                   obj.τs1s2[mask],
                   obj.τs1s3[mask],
                   obj.τs2s3[mask],
                   obj.τm1m2[mask],
                   obj.τm1m3[mask],
                   obj.τm2m3[mask],
                   obj.Δωf[mask]))

end

"""
    readbrainweb(numcomp[, path[, slice]])

This function converts phantom_1.0mm_normal_crisp.rawb (downloaded from
http://brainweb.bic.mni.mcgill.ca/cgi/
brainweb1?alias=phantom_1.0mm_normal_crisp&download=1) into a BrainWeb3Comp
object.

# Input
- `numcomp::Integer`: Number of tissue compartments
- `path::String`: Path to phantom_1.0mm_normal_crisp.rawb
- `slice::Int = 81`: Slice to return (because doing the whole dataset uses a
    lot of memory)
"""
function readbrainweb(numcomp::Integer,
    path::String = modulepath("estimation/data/BrainWeb"),
    slice::Integer = 81
)

    file = joinpath(path, "phantom_1.0mm_normal_crisp.rawb")

    # Make sure BrainWeb data exists
    isfile(file) || throw(ErrorException("Could not find file " *
        "phantom_1.0mm_normal_crisp.rawb at the given path: $path. " *
        "See the README for instructions on how to download the file " *
        "from BrainWeb."))

    # Read the file
    raw = read(file)

    # Reorder the data into a 3D structure
    data = reshape(raw, 181, 217, 181) # [nx,ny,nz]

    # Assign tissue parameters
    if numcomp == 3
        obj = assignvalues3comp(data[:,:,slice])
    elseif numcomp == 300
        obj = assignvalues3compNE(data[:,:,slice])
    elseif numcomp == 4
        obj = assignvalues4comp(data[:,:,slice])
    elseif numcomp == 9
        obj = assignvalues9comp(data[:,:,slice])
    else
        throw(ArgumentError("invalid number of tissue compartments"))
    end

    return obj

end

"""
    assignvalues3comp(id)
"""
function assignvalues3comp(id::Array{UInt8,N}) where {N}

    obj = BrainWeb3Comp(id)

    for index in eachindex(id)
        i = id[index]
        if i == 0 # Background
            obj.tissue[index] = :background
            obj.M0[index]  = 0
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 0
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 0
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 1 # CSF
            obj.tissue[index] = :csf
            obj.M0[index]  = 1
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 3000
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 350
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 2 # Gray Matter
            obj.tissue[index] = :gm
            obj.M0[index]  = 0.86
            obj.ff[index]  = 0.03
            obj.fm[index]  = 0.03
            obj.T1f[index] = 500 # Based on Table 5 in https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.21704
            obj.T1s[index] = 1331
            obj.T1m[index] = 1000
            obj.T2f[index] = 20
            obj.T2s[index] = 80
            obj.T2m[index] = 0.02
            obj.τfs[index] = 20 # 1/5 times WM value because ff is 1/5
            obj.τfm[index] = 10 # Ditto
            obj.Δωf[index] = 5 * 2π
        elseif i == 3 # White Matter
            obj.tissue[index] = :wm
            obj.M0[index]  = 0.77 # From BrainWeb
            obj.ff[index]  = 0.15
            obj.fm[index]  = 0.1
            obj.T1f[index] = 400
            obj.T1s[index] = 832
            obj.T1m[index] = 1000
            obj.T2f[index] = 20
            obj.T2s[index] = 80
            obj.T2m[index] = 0.02
            obj.τfs[index] = 100
            obj.τfm[index] = 50
            obj.Δωf[index] = 15 * 2π
        elseif i == 4 # Fat
            obj.tissue[index] = :fat
            obj.M0[index]  = 1
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 495
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 99
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 5 # Muscle/Skin
            obj.tissue[index] = :muscleskin
            obj.M0[index]  = 1
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 1412
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 50
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 6 # Skin
            obj.tissue[index] = :skin
            obj.M0[index]  = 1
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 3000
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 350
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 7 # Skull
            obj.tissue[index] = :skull
            obj.M0[index]  = 0
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 0
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 0
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 8 # Glial Matter
            obj.tissue[index] = :glial
            obj.M0[index]  = 0.86
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 1331
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 110
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 9 # Meat
            obj.tissue[index] = :meat
            obj.M0[index]  = 0.77
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 832
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 79.6
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 10 # MS Lesion
            obj.tissue[index] = :lesion
            obj.M0[index]  = 0.76
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 1063
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 335
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        else
            throw(ArgumentError("Undefined tissue: id = $i"))
        end
    end

    return obj

end

function assignvalues3compNE(id::Array{UInt8,N}) where {N}

    obj = BrainWeb3Comp(id)

    for index in eachindex(id)
        i = id[index]
        if i == 0 # Background
            obj.tissue[index] = :background
            obj.M0[index]  = 0
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 0
            obj.T1m[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 0
            obj.T2m[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.Δωf[index] = 0
        elseif i == 2 # Gray Matter
            obj.tissue[index] = :gm
            obj.M0[index]  = 0.86
            obj.ff[index]  = 0.03
            obj.fm[index]  = 0.03
            obj.T1f[index] = 500
            obj.T1s[index] = 1331
            obj.T1m[index] = 1000
            obj.T2f[index] = 20
            obj.T2s[index] = 80
            obj.T2m[index] = 0.02
            obj.τfs[index] = Inf
            obj.τfm[index] = Inf
            obj.Δωf[index] = 5 * 2π
        elseif i == 3 # White Matter
            obj.tissue[index] = :wm
            obj.M0[index]  = 0.77
            obj.ff[index]  = 0.15
            obj.fm[index]  = 0.1
            obj.T1f[index] = 400
            obj.T1s[index] = 832
            obj.T1m[index] = 1000
            obj.T2f[index] = 20
            obj.T2s[index] = 80
            obj.T2m[index] = 0.02
            obj.τfs[index] = Inf
            obj.τfm[index] = Inf
            obj.Δωf[index] = 15 * 2π
        else
            obj.tissue[index] = :dontcare
        end
    end

    return obj

end

function assignvalues4comp(id::Array{UInt8,N}) where {N}

    obj = BrainWeb4Comp(id)

    for index in eachindex(id)
        i = id[index]
        if i == 0 # Background
            obj.tissue[index] = :background
            obj.M0[index]  = 0
            obj.ff[index]  = 0
            obj.fm[index]  = 0
            obj.fa[index]  = 0
            obj.T1f[index] = 0
            obj.T1s[index] = 0
            obj.T1m[index] = 0
            obj.T1a[index] = 0
            obj.T2f[index] = 0
            obj.T2s[index] = 0
            obj.T2m[index] = 0
            obj.T2a[index] = 0
            obj.τfs[index] = 0
            obj.τfm[index] = 0
            obj.τfa[index] = 0
            obj.Δωf[index] = 0
        elseif i == 2 # Gray Matter
            obj.tissue[index] = :gm
            obj.M0[index]  = 0.86
            obj.ff[index]  = 0.03
            obj.fm[index]  = 0.03
            obj.fa[index]  = 0.5
            obj.T1f[index] = 500
            obj.T1s[index] = 1331
            obj.T1m[index] = 1000
            obj.T1a[index] = 3000
            obj.T2f[index] = 20
            obj.T2s[index] = 80
            obj.T2m[index] = 0.02
            obj.T2a[index] = 150
            obj.τfs[index] = 40
            obj.τfm[index] = 10
            obj.τfa[index] = 40
            obj.Δωf[index] = 5 * 2π
        elseif i == 3 # White Matter
            obj.tissue[index] = :wm
            obj.M0[index]  = 0.77
            obj.ff[index]  = 0.15
            obj.fm[index]  = 0.1
            obj.fa[index]  = 0.3
            obj.T1f[index] = 400
            obj.T1s[index] = 832
            obj.T1m[index] = 1000
            obj.T1a[index] = 3000
            obj.T2f[index] = 20
            obj.T2s[index] = 80
            obj.T2m[index] = 0.02
            obj.T2a[index] = 150
            obj.τfs[index] = 200
            obj.τfm[index] = 50
            obj.τfa[index] = 200
            obj.Δωf[index] = 15 * 2π
        else
            obj.tissue[index] = :dontcare
        end
    end

    return obj

end

function assignvalues9comp(id::Array{UInt8,N}) where {N}

    obj = BrainWeb9Comp(id)

    for index in eachindex(id)
        i = id[index]
        if i == 0 # Background
            obj.tissue[index] = :background
            obj.M0[index]     = 0
            obj.ff1[index]    = 0
            obj.ff2[index]    = 0
            obj.ff3[index]    = 0
            obj.fs1[index]    = 0
            obj.fs2[index]    = 0
            obj.fs3[index]    = 0
            obj.fm1[index]    = 0
            obj.fm2[index]    = 0
            obj.fm3[index]    = 0
            obj.T1f1[index]   = 0
            obj.T1f2[index]   = 0
            obj.T1f3[index]   = 0
            obj.T1s1[index]   = 0
            obj.T1s2[index]   = 0
            obj.T1s3[index]   = 0
            obj.T1m1[index]   = 0
            obj.T1m2[index]   = 0
            obj.T1m3[index]   = 0
            obj.T2f1[index]   = 0
            obj.T2f2[index]   = 0
            obj.T2f3[index]   = 0
            obj.T2s1[index]   = 0
            obj.T2s2[index]   = 0
            obj.T2s3[index]   = 0
            obj.T2m1[index]   = 0
            obj.T2m2[index]   = 0
            obj.T2m3[index]   = 0
            obj.τf1f2[index]  = 0
            obj.τf1f3[index]  = 0
            obj.τf1s1[index]  = 0
            obj.τf1s2[index]  = 0
            obj.τf1s3[index]  = 0
            obj.τf1m1[index]  = 0
            obj.τf1m2[index]  = 0
            obj.τf1m3[index]  = 0
            obj.τf2f3[index]  = 0
            obj.τf2s1[index]  = 0
            obj.τf2s2[index]  = 0
            obj.τf2s3[index]  = 0
            obj.τf2m1[index]  = 0
            obj.τf2m2[index]  = 0
            obj.τf2m3[index]  = 0
            obj.τf3s1[index]  = 0
            obj.τf3s2[index]  = 0
            obj.τf3s3[index]  = 0
            obj.τf3m1[index]  = 0
            obj.τf3m2[index]  = 0
            obj.τf3m3[index]  = 0
            obj.τs1s2[index]  = 0
            obj.τs1s3[index]  = 0
            obj.τs2s3[index]  = 0
            obj.τm1m2[index]  = 0
            obj.τm1m3[index]  = 0
            obj.τm2m3[index]  = 0
            obj.Δωf[index]    = 0
        elseif i == 2 # Gray Matter
            obj.tissue[index] = :gm
            obj.M0[index]     = 0.86
            obj.ff1[index]    = 0.0075
            obj.ff2[index]    = 0.015
            obj.ff3[index]    = 0.0075
            obj.fs1[index]    = 0.235
            obj.fs2[index]    = 0.47
            obj.fs3[index]    = 0.235
            obj.fm1[index]    = 0.0075
            obj.fm2[index]    = 0.015
            obj.fm3[index]    = 0.0075
            obj.T1f1[index]   = 500
            obj.T1f2[index]   = 500
            obj.T1f3[index]   = 500
            obj.T1s1[index]   = 1331
            obj.T1s2[index]   = 1331
            obj.T1s3[index]   = 1331
            obj.T1m1[index]   = 1000
            obj.T1m2[index]   = 1000
            obj.T1m3[index]   = 1000
            obj.T2f1[index]   = 16
            obj.T2f2[index]   = 20
            obj.T2f3[index]   = 24
            obj.T2s1[index]   = 64
            obj.T2s2[index]   = 80
            obj.T2s3[index]   = 96
            obj.T2m1[index]   = 0.01
            obj.T2m2[index]   = 0.02
            obj.T2m3[index]   = 0.03
            obj.τf1f2[index]  = 10
            obj.τf1f3[index]  = 10
            obj.τf1s1[index]  = 40
            obj.τf1s2[index]  = 80
            obj.τf1s3[index]  = 80
            obj.τf1m1[index]  = 20
            obj.τf1m2[index]  = 40
            obj.τf1m3[index]  = 40
            obj.τf2f3[index]  = 10
            obj.τf2s1[index]  = 80
            obj.τf2s2[index]  = 40
            obj.τf2s3[index]  = 80
            obj.τf2m1[index]  = 40
            obj.τf2m2[index]  = 20
            obj.τf2m3[index]  = 40
            obj.τf3s1[index]  = 80
            obj.τf3s2[index]  = 80
            obj.τf3s3[index]  = 40
            obj.τf3m1[index]  = 40
            obj.τf3m2[index]  = 40
            obj.τf3m3[index]  = 20
            obj.τs1s2[index]  = 10
            obj.τs1s3[index]  = 10
            obj.τs2s3[index]  = 10
            obj.τm1m2[index]  = 10
            obj.τm1m3[index]  = 10
            obj.τm2m3[index]  = 10
            obj.Δωf[index]    = 5 * 2π
        elseif i == 3 # White Matter
            obj.tissue[index] = :wm
            obj.M0[index]     = 0.77
            obj.ff1[index]    = 0.0375
            obj.ff2[index]    = 0.075
            obj.ff3[index]    = 0.0375
            obj.fs1[index]    = 0.1875
            obj.fs2[index]    = 0.375
            obj.fs3[index]    = 0.1875
            obj.fm1[index]    = 0.025
            obj.fm2[index]    = 0.05
            obj.fm3[index]    = 0.025
            obj.T1f1[index]   = 400
            obj.T1f2[index]   = 400
            obj.T1f3[index]   = 400
            obj.T1s1[index]   = 832
            obj.T1s2[index]   = 832
            obj.T1s3[index]   = 832
            obj.T1m1[index]   = 1000
            obj.T1m2[index]   = 1000
            obj.T1m3[index]   = 1000
            obj.T2f1[index]   = 16
            obj.T2f2[index]   = 20
            obj.T2f3[index]   = 24
            obj.T2s1[index]   = 64
            obj.T2s2[index]   = 80
            obj.T2s3[index]   = 96
            obj.T2m1[index]   = 0.01
            obj.T2m2[index]   = 0.02
            obj.T2m3[index]   = 0.03
            obj.τf1f2[index]  = 10
            obj.τf1f3[index]  = 10
            obj.τf1s1[index]  = 200
            obj.τf1s2[index]  = 400
            obj.τf1s3[index]  = 400
            obj.τf1m1[index]  = 100
            obj.τf1m2[index]  = 200
            obj.τf1m3[index]  = 200
            obj.τf2f3[index]  = 10
            obj.τf2s1[index]  = 400
            obj.τf2s2[index]  = 200
            obj.τf2s3[index]  = 400
            obj.τf2m1[index]  = 200
            obj.τf2m2[index]  = 100
            obj.τf2m3[index]  = 200
            obj.τf3s1[index]  = 400
            obj.τf3s2[index]  = 400
            obj.τf3s3[index]  = 200
            obj.τf3m1[index]  = 200
            obj.τf3m2[index]  = 200
            obj.τf3m3[index]  = 100
            obj.τs1s2[index]  = 10
            obj.τs1s3[index]  = 10
            obj.τs2s3[index]  = 10
            obj.τm1m2[index]  = 10
            obj.τm1m3[index]  = 10
            obj.τm2m3[index]  = 10
            obj.Δωf[index]    = 15 * 2π
        else
            obj.tissue[index] = :dontcare
        end
    end

    return obj

end
