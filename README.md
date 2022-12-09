# BioFibers.jl

This module provides two core types for building data structures that store physical properties pertinent to  modeling biologic fibers.

To use this module you will need to add the following Julia packages to your implementation of `BioFibers`:

```
using Pkg
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add(url = "https://github.com/AlanFreed/MutableTypes.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalSystemsOfUnits.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalScalars.jl")
```

## BioFiber

Type `BioFiber` holds the physical data that describes a generic biologic fiber.

```
struct BioFiber
    # structural property
    ruptured::MBool       # specifies if fiber has ruptured, or not

    # physical property
    ϵf::PhysicalScalar    # fracture strain
    ρ::PhysicalScalar     # mass density                     g/cm³

    # thermal properties
    α::PhysicalScalar     # thermal strain coefficient      1/°C
    cₚ::PhysicalScalar    # specific heat at constant P      erg/(g⋅°C) = cm²/(s²⋅°C)

    # Fiber properties in a reference configuration κᵣ:

    # dimensional properties
    Lᵣ::PhysicalScalar    # reference length                 cm
    Aᵣ::PhysicalScalar    # reference cross-sectional area   cm²

    # thermodynamic conjugate fields: causes
    θᵣ::PhysicalScalar    # reference temperature            °C
    σᵣ::PhysicalScalar    # reference (residual) stress      g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    ηᵣ::PhysicalScalar    # reference entropy per unit mass  cm²/(s²⋅°C)
    ϵᵣ::PhysicalScalar    # reference strain

    # Fiber properties in the initial configurations κ₀:

    # dimensional properties
    L₀::PhysicalScalar    # initial length                   cm
    A₀::PhysicalScalar    # initial cross-sectional area     cm²

    # thermodynamic conjugate fields: causes
    θ₀::PhysicalScalar    # initial temperature              °C
    σ₀::PhysicalScalar    # initial stress                   g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    η₀::PhysicalScalar    # initial entropy per unit mass    cm²/(s²⋅°C)
    ϵ₀::PhysicalScalar    # initial strain
end
```

### Methods

To make a shallow copy of an instance of type `BioFiber`:
```
function Base.:(copy)(f::BioFiber)::BioFiber
```

To make a deep copy of an instance of type `BioFiber`:
```
function Base.:(deepcopy)(f::BioFiber)::BioFiber
```

To convert an instance of type `BioFiber` into a `String`:
```
function toString(f::BioFiber;
                  format::Char='E',
                  precision::Int=5,
                  aligned::Bool=false)::String
```
where the keyword `format` is a character that, whenever its value is 'E' or 'e', will represent scalar fields in scientific notation; otherwise, they will be represented in fixed-point notation. Keyword `precision` specifies the number of significant digits to be shown, which can accept values from the set \{3…7\}. Keyword `aligned`, when set to `true`, will add a white space in front of any non-negative scalar string representation, e.g., this could be useful when printing out a matrix of scalars; otherwise, there is no leading white space in its string representation, which is the default.


## AlveolarChord

Alveolar chords are comprised of collagen and elastin fibers loaded in parallel with weak mechanical coupling between them. Type `AlveolarChord` holds the physical data that describes a generic alveolar chord.

```
struct AlveolarChord
    # fibers of the chord
    fᶜ::BioFiber          # collagen fiber in the alveolar chord
    fᵉ::BioFiber          # elastin  fiber in the alveolar chord

    # structural property
    ruptured::MBool       # specifies if fiber has ruptured, or not

    # physical properties
    ρ::PhysicalScalar     # mass density                     g/cm³
    ϕᶜ::PhysicalScalar    # volume fraction of collagen

    # Chordal properties in a reference configuration κᵣ:

    # dimensional properties
    Lᵣ::PhysicalScalar    # reference length                 cm
    Aᵣ::PhysicalScalar    # reference cross-sectional area   cm²

    # thermodynamic conjugate fields: causes
    θᵣ::PhysicalScalar    # reference temperature            °C
    σᵣ::PhysicalScalar    # reference (residual) stress      g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    ηᵣ::PhysicalScalar    # reference entropy per unit mass  cm²/(s²⋅°C)
    ϵᵣ::PhysicalScalar    # reference strain

    # Chordal properties in the initial configurations κ₀:

    # dimensional properties
    L₀::PhysicalScalar    # initial length                   cm
    A₀::PhysicalScalar    # initial cross-sectional area     cm²

    # thermodynamic conjugate fields: causes
    θ₀::PhysicalScalar    # initial temperature              °C
    σ₀::PhysicalScalar    # initial stress                   g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    η₀::PhysicalScalar    # initial entropy per unit mass    cm²/(s²⋅°C)
    ϵ₀::PhysicalScalar    # initial strain
end
```

### Methods

To make a shallow copy of an instance of type `AlveolarChord`:
```
function Base.:(copy)(ac::AlveolarChord)::AlveolarChord
```

To make a deep copy of an instance of type `AlveolarChord`:
```
function Base.:(deepcopy)(ac::AlveolarChord)::AlveolarChord
```

To convert an instance of type `AlveolarChord` into a `String`:
```
function toString(ac::AlveolarChord;
                  format::Char='E',
                  precision::Int=5,
                  aligned::Bool=true)::String
```
where the keyword `format` is a character that, whenever its value is 'E' or 'e', will represents scalars in scientific notation; otherwise, they will be represented in fixed-point notation. Keyword `precision` specifies the number of significant digits to be shown, which can accept values from the set \{3…7\}. Keyword `aligned`, when set to `true`, will add a white space in front of any non-negative scalar string representation, e.g., this could be useful when printing out a matrix of scalars; otherwise, there is no leading white space in its string representation, which is the default.

### Consturctors

To create a new instance of type `BioFiber` for a collagen fiber in an alveolar chord:
```
function newCollagenFiber(Lᵣ::PhysicalScalar, L₀::PhysicalScalar)::BioFiber
```
where `Lᵣ` is its length in the reference (strain-free) configuration κᵣ, while `L₀` is its length in an initial configuration κ₀ of analysis, which is typically distinct from κᵣ.

To create a new instance of type `BioFiber` for an elastin fiber in an alveolar chord:
```
function newElastinFiber(Lᵣ::PhysicalScalar, L₀::PhysicalScalar)::BioFiber
```
where `Lᵣ` is its length in the reference (strain-free) configuration κᵣ, while `L₀` is its length in an initial configuration κ₀ of analysis, which is typically distinct from κᵣ.

To create a new instance of type `AlveolarChord`:
```
function newAlveolarChord(Lᵣ::PhysicalScalar, L₀::PhysicalScalar)::AlveolarChord
```
where `Lᵣ` is its chordal length in the reference (strain-free) configuration κᵣ, while `L₀` is its chordal length in an initial configuration κ₀ of analysis, which is typically distinct from κᵣ.