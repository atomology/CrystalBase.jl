using PeriodicTable: elements

export atomic_number, atomic_symbol

"""
    $(SIGNATURES)

Get atom number from symbol.

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
    return [findfirst(x -> x == string(s), table) for s in symbols]
end

"""
    $(SIGNATURES)

Get atomic symbol from number.

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
    return [elements[n].symbol for n in numbers]
end
