"""
    getdata(doi, filename)

Download the data from the given DOI. This function programmatically "clicks" on
the download button found at the Deep Blue Data website given by `doi`.

# Arguments
- `doi::AbstractString`: DOI from which to download the data
- `filename::AbstractString`: Name (and path) for the downloaded data
"""
function getdata(doi::AbstractString, filename::AbstractString)

    response = HTTP.get(url)
    doc = parsehtml(String(response.body))
    links = Vector{String}()
    findlinks!(links, doc.root)
    isempty(links) && throw(ErrorException("No download links were found at " *
                                           "the given DOI: $doi"))
    length(links) > 1 && throw(ErrorExecption("Multiple download links found " *
                                              "at the given DOI: $doi"))
    link = joinpath("https://deepblue.lib.umich.edu", links[])
    print("Downloading data from $doi to $filename...")
    download(link, filename)
    println("done!")

end

"""
    findlinks!(links, elem)

Recursively search `elem` and all its children for download links, adding them
to `links`.
"""
function findlinks!(links::AbstractVector{<:AbstractString}, elem::HTMLElement)

    if haskey(attrs(elem), "id") && getattr(elem, "id") == "file_download"
        haskey(attrs(elem), "href") && push!(links, getattr(elem, "href"))
    end
    c = children(elem)
    foreach(el -> findlinks!(links, el),
            (c[i] for i = 1:length(c) if c[i] isa HTMLElement))

end
