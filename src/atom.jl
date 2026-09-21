using PeriodicTable: elements

export atomic_number, atomic_symbol

"""
    $(SIGNATURES)

Get atom number from symbol. The dummy symbol `"X"` maps to `0`, see
[`DUMMY_SYMBOL`](@ref).

# Examples
```jldoctest atomic_number; setup = :(using CrystalBase)
atomic_number("Si")
# output
14
```

```jldoctest atomic_number
atomic_number(:O)
# output
8
```

```jldoctest atomic_number
atomic_number(["Si", :O])
# output
2-element Vector{Int64}:
 14
  8
```
"""
function atomic_number(symbol::Union{AbstractString, Symbol})
    return atomic_number([symbol])[1]
end

function atomic_number(symbols::AbstractVector)
    table = [e.symbol for e in elements]
    return map(symbols) do s
        string(s) == DUMMY_SYMBOL && return 0
        findfirst(x -> x == string(s), table)
    end
end

"""
Symbol of the dummy element with atomic number `0`, for sites of synthetic
models that carry no chemical identity.
"""
const DUMMY_SYMBOL = "X"

"""
    $(SIGNATURES)

Get atomic symbol from number. `0` maps to the dummy symbol `"X"`.

# Examples
```jldoctest atomic_symbol; setup = :(using CrystalBase)
atomic_symbol(14)
# output
"Si"
```

```jldoctest atomic_symbol
atomic_symbol([14, 8])
# output
2-element Vector{String}:
 "Si"
 "O"
```
"""
function atomic_symbol(number::Integer)
    return atomic_symbol([number])[1]
end

function atomic_symbol(numbers::AbstractVector)
    return [n == 0 ? DUMMY_SYMBOL : elements[n].symbol for n in numbers]
end

export label_symbol, formula

"""
    $(SIGNATURES)

Element symbol at the head of a species label, e.g. `"Fe1"` → `"Fe"`.

A species label starts with an element symbol optionally followed by a
disambiguating suffix (`Fe1`, `Fe_up`, `Fe3+`). The longest valid prefix wins,
so two-letter symbols are preferred over one (`"Co"` is cobalt, not carbon).
The dummy symbol `"X"` (atomic number 0) is accepted, so `"X1"` labels a
site of a synthetic model.

Throws `ArgumentError` if the label does not begin with a known element.

# Examples
```jldoctest label_symbol; setup = :(using CrystalBase)
label_symbol("Fe1")
# output
"Fe"
```

```jldoctest label_symbol
label_symbol.(["Si", "O2", "Cl_dn"])
# output
3-element Vector{String}:
 "Si"
 "O"
 "Cl"
```
"""
function label_symbol(label::AbstractString)
    s = String(label)
    for width in (2, 1)
        ncodeunits(s) >= width || continue
        head = first(s, width)
        isnothing(atomic_number(head)) || return head
    end
    throw(ArgumentError("label \"$label\" does not begin with a known element symbol"))
end

"""
    $(SIGNATURES)

Render element `symbols` (one entry per atom) as a chemical formula string.

A count of 1 is written bare (`SiO2`, not `Si1O2`). With `reduced = true` the
counts are divided by their greatest common divisor first (`Si2O4` → `SiO2`).

`order` selects the element ordering:
- `:hill` (default): Hill order, the CAS/database convention. With carbon
  present, `C` then `H` then the rest alphabetical; without carbon, everything
  alphabetical (`O2Si`).
- `:alpha`: plain alphabetical.

# Examples
```jldoctest formula; setup = :(using CrystalBase)
formula(["Si", "O", "O"])
# output
"O2Si"
```

```jldoctest formula
formula(["C", "O", "H", "H", "H", "H"])
# output
"CH4O"
```

```jldoctest formula
formula(["Fe", "Fe", "Fe", "Fe", "O", "O", "O", "O", "O", "O"]; reduced = true)
# output
"Fe2O3"
```
"""
function formula(symbols::AbstractVector{<:AbstractString}; order::Symbol = :hill, reduced::Bool = false)
    counts = OrderedDict{String, Int}()
    for s in symbols
        counts[String(s)] = get(counts, String(s), 0) + 1
    end
    if reduced && !isempty(counts)
        g = reduce(gcd, values(counts))
        for (k, v) in counts
            counts[k] = v ÷ g
        end
    end
    ordered = _formula_order(collect(keys(counts)), order)
    return join(string(e, counts[e] > 1 ? counts[e] : "") for e in ordered)
end

function _formula_order(elements::Vector{String}, order::Symbol)
    if order === :alpha
        return sort(elements)
    elseif order === :hill
        if "C" in elements
            head = ["C"]
            "H" in elements && push!(head, "H")
            return vcat(head, sort(filter(e -> e ∉ ("C", "H"), elements)))
        end
        return sort(elements)
    end
    throw(ArgumentError("unknown formula order: $order (expected :hill or :alpha)"))
end
