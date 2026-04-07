using PeriodicTable: elements

export atomic_number, atomic_symbol

"""
    $(SIGNATURES)

Get atom number from symbol.
"""
function atomic_number(symbol::Union{AbstractString, Symbol})
    return atomic_number([symbol])[1]
end

"""
    $(SIGNATURES)

Get atom number from symbol.
"""
function atomic_number(symbols::AbstractVector)
    table = [e.symbol for e in elements]
    return [findfirst(x -> x == string(s), table) for s in symbols]
end

"""
    $(SIGNATURES)

Get atomic symbol from number.
"""
function atomic_symbol(number::Integer)
    return atomic_symbol([number])[1]
end

"""
    $(SIGNATURES)
Get atomic symbol from number.
"""
function atomic_symbol(numbers::AbstractVector)
    return [elements[n].symbol for n in numbers]
end
