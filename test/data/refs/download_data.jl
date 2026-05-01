# Pulls each survival dataset from the upstream Rdatasets repository
# (https://github.com/vincentarelbundock/Rdatasets) into `<name>_data.csv`,
# which both `generate_refs.R` and `runtests.jl` consume. Run after the
# upstream repo is updated; commit the resulting CSVs.
#
# Usage:
#   julia test/data/refs/download_data.jl
#
# Note: `csv/survival/lung.csv` and `csv/survival/rats.csv` in upstream are
# clobbered by same-named files from other R packages (lung→data column only,
# rats→aod's diet experiment). We pull lung from `cancer.csv` (identical to
# survival's lung) and skip rats entirely.

using Downloads

const BASE  = "https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/survival/"
const HERE  = @__DIR__

# (output name, upstream filename)
const FILES = [
    "lung"      => "cancer.csv",
    "leukemia"  => "leukemia.csv",
    "mgus"      => "mgus.csv",
    "nwtco"     => "nwtco.csv",
    "ovarian"   => "ovarian.csv",
    "pbc"       => "pbc.csv",
    "stanford2" => "stanford2.csv",
    "veteran"   => "veteran.csv",
    "kidney"    => "kidney.csv",
    "colon"     => "colon.csv",
]

# Drop the `rownames` index column and rewrite `.` → `_` in the header so the
# columns are valid Julia identifiers usable in `@formula`.
function clean_csv!(path)
    lines = readlines(path)
    isempty(lines) && error("$(path) is empty")
    header = split(lines[1], ',')
    drop = findall(==("rownames"), header)
    keep = setdiff(eachindex(header), drop)
    new_header = replace.(header[keep], '.' => '_')
    open(path, "w") do io
        println(io, join(new_header, ','))
        for ln in @view lines[2:end]
            fields = split(ln, ',')
            println(io, join(fields[keep], ','))
        end
    end
end

for (name, src) in FILES
    out = joinpath(HERE, "$(name)_data.csv")
    Downloads.download(BASE * src, out)
    clean_csv!(out)
    println("wrote ", relpath(out))
end
