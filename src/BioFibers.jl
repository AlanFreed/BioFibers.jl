module BioFibers

using
    Distributions,
    Random

import
    MutableTypes:
        MBool,
        get                 as Mget

import
    PhysicalSystemsOfUnits:
        CGS,
        BARYE,
        CENTIGRADE,
        CENTIMETER,
        CGS_AREA            as AREA,
        CGS_DIMENSIONLESS   as DIMENSIONLESS,
        CGS_ENTROPYperMASS  as ENTROPYperUnitMASS,
        CGS_MASS_DENSITY    as MASSDENSITY,
        CGS_MODULUS         as MODULUS

import
    PhysicalFields:
        PhysicalScalar      as PhyScalar,
        newPhysicalScalar   as newScalar,
        StoString

import
    PhysicalScalars:
        ==, ≈, ≠, <, ≤, ≥, >, +, -, *, /, ^,
        get                 as Sget,
        set!,
        exp,
        log

export
    # types
    BioFiber,
    AlveolarChord,
    # methods
    copy,
    deepcopy,
    toString,
    # function
    newAlveolarChord

const THERMALSTRAIN = CGS(0, 0, 0, -1)

struct BioFiber
    # structural property
    ruptured::MBool  # specifies if fiber has ruptured, or not

    # physical property
    ϵf::PhyScalar    # fracture strain
    ρ::PhyScalar     # mass density                     g/cm³

    # thermal properties
    α::PhyScalar     # thermal strain coefficient      1/°C
    cₚ::PhyScalar    # specific heat at constant P      erg/(g⋅°C) = cm²/(s²⋅°C)

    # Fiber properties in a reference configuration κᵣ:

    # dimensional properties
    Lᵣ::PhyScalar    # reference length                 cm
    Aᵣ::PhyScalar    # reference cross-sectional area   cm²

    # thermodynamic conjugate fields: causes
    θᵣ::PhyScalar    # reference temperature            °C
    σᵣ::PhyScalar    # reference (residual) stress      g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    ηᵣ::PhyScalar    # reference entropy per unit mass  cm²/(s²⋅°C)
    ϵᵣ::PhyScalar    # reference strain

    # Fiber properties in the initial configurations κ₀:

    # dimensional properties
    L₀::PhyScalar    # initial length                   cm
    A₀::PhyScalar    # initial cross-sectional area     cm²

    # thermodynamic conjugate fields: causes
    θ₀::PhyScalar    # initial temperature              °C
    σ₀::PhyScalar    # initial stress                   g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    η₀::PhyScalar    # initial entropy per unit mass    cm²/(s²⋅°C)
    ϵ₀::PhyScalar    # initial strain
end # BioFiber

struct AlveolarChord
    # fibers of the chord
    fᶜ::BioFiber     # collagen fiber in the alveolar chord
    fᵉ::BioFiber     # elastin  fiber in the alveolar chord

    # structural property
    ruptured::MBool  # specifies if fiber has ruptured, or not

    # physical properties
    ρ::PhyScalar     # mass density                     g/cm³
    ϕᶜ::PhyScalar    # volume fraction of collagen

    # Chordal properties in a reference configuration κᵣ:

    # dimensional properties
    Lᵣ::PhyScalar    # reference length                 cm
    Aᵣ::PhyScalar    # reference cross-sectional area   cm²

    # thermodynamic conjugate fields: causes
    θᵣ::PhyScalar    # reference temperature            °C
    σᵣ::PhyScalar    # reference (residual) stress      g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    ηᵣ::PhyScalar    # reference entropy per unit mass  cm²/(s²⋅°C)
    ϵᵣ::PhyScalar    # reference strain

    # Chordal properties in the initial configurations κ₀:

    # dimensional properties
    L₀::PhyScalar    # initial length                   cm
    A₀::PhyScalar    # initial cross-sectional area     cm²

    # thermodynamic conjugate fields: causes
    θ₀::PhyScalar    # initial temperature              °C
    σ₀::PhyScalar    # initial stress                   g/(cm⋅s²)

    # thermodynamic conjugate fields: effects
    η₀::PhyScalar    # initial entropy per unit mass    cm²/(s²⋅°C)
    ϵ₀::PhyScalar    # initial strain
end # AlveolarChord

# methods

function Base.:(copy)(f::BioFiber)::BioFiber
    ruptured = copy(f.ruptured)
    ϵf = copy(f.ϵf)
    ρ  = copy(f.ρ)
    α  = copy(f.α)
    cₚ = copy(f.cₚ)
    Lᵣ = copy(f.Lᵣ)
    Aᵣ = copy(f.Aᵣ)
    θᵣ = copy(f.θᵣ)
    σᵣ = copy(f.σᵣ)
    ηᵣ = copy(f.ηᵣ)
    ϵᵣ = copy(f.ϵᵣ)
    L₀ = copy(f.L₀)
    A₀ = copy(f.A₀)
    θ₀ = copy(f.θ₀)
    σ₀ = copy(f.σ₀)
    η₀ = copy(f.η₀)
    ϵ₀ = copy(f.ϵ₀)
    return BioFiber(ruptured, ϵf, ρ, α, cₚ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end

function Base.:(deepcopy)(f::BioFiber)::BioFiber
    ruptured = deepcopy(f.ruptured)
    ϵf = deepcopy(f.ϵf)
    ρ  = deepcopy(f.ρ)
    α  = deepcopy(f.α)
    cₚ = deepcopy(f.cₚ)
    Lᵣ = deepcopy(f.Lᵣ)
    Aᵣ = deepcopy(f.Aᵣ)
    θᵣ = deepcopy(f.θᵣ)
    σᵣ = deepcopy(f.σᵣ)
    ηᵣ = deepcopy(f.ηᵣ)
    ϵᵣ = deepcopy(f.ϵᵣ)
    L₀ = deepcopy(f.L₀)
    A₀ = deepcopy(f.A₀)
    θ₀ = deepcopy(f.θ₀)
    σ₀ = deepcopy(f.σ₀)
    η₀ = deepcopy(f.η₀)
    ϵ₀ = deepcopy(f.ϵ₀)
    return BioFiber(ruptured, ϵf, ρ, α, cₚ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end

function Base.:(copy)(ac::AlveolarChord)::AlveolarChord
    fᶜ = copy(ac.fᶜ)
    fᵉ = copy(ac.fᵉ)
    ruptured = copy(ac.ruptured)
    ρ  = copy(ac.ρ)
    ϕᶜ = copy(ac.ϕᶜ)
    Lᵣ = copy(ac.Lᵣ)
    Aᵣ = copy(ac.Aᵣ)
    θᵣ = copy(ac.θᵣ)
    σᵣ = copy(ac.σᵣ)
    ηᵣ = copy(ac.ηᵣ)
    ϵᵣ = copy(ac.ϵᵣ)
    L₀ = copy(ac.L₀)
    A₀ = copy(ac.A₀)
    θ₀ = copy(ac.θ₀)
    σ₀ = copy(ac.σ₀)
    η₀ = copy(ac.η₀)
    ϵ₀ = copy(ac.ϵ₀)
    return AlveolarChord(fᶜ, fᵉ, ruptured, ρ, ϕᶜ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end

function Base.:(deepcopy)(ac::AlveolarChord)::AlveolarChord
    fᶜ = deepcopy(ac.fᶜ)
    fᵉ = deepcopy(ac.fᵉ)
    ruptured = deepcopy(ac.ruptured)
    ρ  = deepcopy(ac.ρ)
    ϕᶜ = deepcopy(ac.ϕᶜ)
    Lᵣ = deepcopy(ac.Lᵣ)
    Aᵣ = deepcopy(ac.Aᵣ)
    θᵣ = deepcopy(ac.θᵣ)
    σᵣ = deepcopy(ac.σᵣ)
    ηᵣ = deepcopy(ac.ηᵣ)
    ϵᵣ = deepcopy(ac.ϵᵣ)
    L₀ = deepcopy(ac.L₀)
    A₀ = deepcopy(ac.A₀)
    θ₀ = deepcopy(ac.θ₀)
    σ₀ = deepcopy(ac.σ₀)
    η₀ = deepcopy(ac.η₀)
    ϵ₀ = deepcopy(ac.ϵ₀)
    return AlveolarChord(fᶜ, fᵉ, ruptured, ρ, ϕᶜ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end

function toString(f::BioFiber;
                  format::Char='E',
                  precision::Int=5,
                  aligned::Bool=false)::String
    if Mget(f.ruptured)
        s = "This fiber has ruptured.\n"
    else
        s = "This fiber has material properties of:\n"
        s = string(s, "   physical:\n")
        s = string(s, "         ϵf = ", StoString(f.ϵf; format, precision, aligned), "\n")
        s = string(s, "         ρ  = ", StoString(f.ρ; format, precision, aligned), "\n")
        s = string(s, "   thermal:\n")
        s = string(s, "         α  = ", StoString(f.α; format, precision, aligned), "\n")
        s = string(s, "         cₚ = ", StoString(f.cₚ; format, precision, aligned), "\n")
        s = string(s, "In its reference configuration κᵣ,\n")
        s = string(s, "   its physical dimensions are:\n")
        s = string(s, "         Lᵣ = ", StoString(f.Lᵣ; format, precision, aligned), "\n")
        s = string(s, "         Aᵣ = ", StoString(f.Aᵣ; format, precision, aligned), "\n")
        s = string(s, "   its thermodynamic conjugate fields are:\n")
        s = string(s, "      causes:\n")
        s = string(s, "         θᵣ = ", StoString(f.θᵣ; format, precision, aligned), "\n")
        s = string(s, "         σᵣ = ", StoString(f.σᵣ; format, precision, aligned), "\n")
        s = string(s, "      effects:\n")
        s = string(s, "         ηᵣ = ", StoString(f.ηᵣ; format, precision, aligned), "\n")
        s = string(s, "         ϵᵣ = ", StoString(f.ϵᵣ; format, precision, aligned), "\n")
        s = string(s, "In its initial configuration κ₀,\n")
        s = string(s, "   its physical dimensions are:\n")
        s = string(s, "         L₀ = ", StoString(f.L₀; format, precision, aligned), "\n")
        s = string(s, "         A₀ = ", StoString(f.A₀; format, precision, aligned), "\n")
        s = string(s, "   its thermodynamic conjugate fields are:\n")
        s = string(s, "      causes:\n")
        s = string(s, "         θ₀ = ", StoString(f.θ₀; format, precision, aligned), "\n")
        s = string(s, "         σ₀ = ", StoString(f.σ₀; format, precision, aligned), "\n")
        s = string(s, "      effects:\n")
        s = string(s, "         η₀ = ", StoString(f.η₀; format, precision, aligned), "\n")
        s = string(s, "         ϵ₀ = ", StoString(f.ϵ₀; format, precision, aligned), "\n")
    end
    return s
end

function toString(ac::AlveolarChord;
                  format::Char='E',
                  precision::Int=5,
                  aligned::Bool=false)::String
    if Mget(ac.ruptured)
        s = "This alveolar chord has ruptured.\n"
    else
        s = "This alveolar chord has a collagen fiber with properties:\n"
        s = string(s, toString(ac.fᶜ; format, precision, aligned))
        s = string(s, "This alveolar chord also has an elastin fiber with properties:\n")
        s = string(s, toString(ac.fᵉ; format, precision, aligned))
        s = string(s, "The alveolar chord has material properties of:\n")
        s = string(s, "   physical:\n")
        s = string(s, "         ϕᶜ = ", StoString(ac.ϕᶜ; format, precision, aligned), "\n")
        s = string(s, "         ρ  = ", StoString(ac.ρ; format, precision, aligned), "\n")
        s = string(s, "In its reference configuration κᵣ,\n")
        s = string(s, "   its physical dimensions are:\n")
        s = string(s, "         Lᵣ = ", StoString(ac.Lᵣ; format, precision, aligned), "\n")
        s = string(s, "         Aᵣ = ", StoString(ac.Aᵣ; format, precision, aligned), "\n")
        s = string(s, "   its thermodynamic conjugate fields are:\n")
        s = string(s, "      causes:\n")
        s = string(s, "         θᵣ = ", StoString(ac.θᵣ; format, precision, aligned), "\n")
        s = string(s, "         σᵣ = ", StoString(ac.σᵣ; format, precision, aligned), "\n")
        s = string(s, "      effects:\n")
        s = string(s, "         ηᵣ = ", StoString(ac.ηᵣ; format, precision, aligned), "\n")
        s = string(s, "         ϵᵣ = ", StoString(ac.ϵᵣ; format, precision, aligned), "\n")
        s = string(s, "In its initial configuration κ₀,\n")
        s = string(s, "   its physical dimensions are:\n")
        s = string(s, "         L₀ = ", StoString(ac.L₀; format, precision, aligned), "\n")
        s = string(s, "         A₀ = ", StoString(ac.A₀; format, precision, aligned), "\n")
        s = string(s, "   its thermodynamic conjugate fields are:\n")
        s = string(s, "      causes:\n")
        s = string(s, "         θ₀ = ", StoString(ac.θ₀; format, precision, aligned), "\n")
        s = string(s, "         σ₀ = ", StoString(ac.σ₀; format, precision, aligned), "\n")
        s = string(s, "      effects:\n")
        s = string(s, "         η₀ = ", StoString(ac.η₀; format, precision, aligned), "\n")
        s = string(s, "         ϵ₀ = ", StoString(ac.ϵ₀; format, precision, aligned), "\n")
    end
    return s
end

# functions

function _collagen(Lᵣ, L₀::PhyScalar)::BioFiber
    # The fiber is originally intact.
    ruptured = MBool(false)

    # Physical properties include a rupture strain of
    μ = 0.25    # mean value for the rupture strain
    σ = 0.025   # standard deviation of the rupture strain
    ϵr = rand(Normal(μ, σ))
    while ϵr < 0.1
        ϵr = rand(Normal(μ, σ))
    end
    ϵf = newScalar(DIMENSIONLESS)
    set!(ϵf, ϵr)
    # and the mass density for hydrated collagen in g/cm³ of
    ρ = newScalar(MASSDENSITY)
    set!(ρ, 1.34)

    # Thermal properties include a thermal strain in 1/°C of
    α = newScalar(THERMALSTRAIN)
    set!(α, 0.001)
    # and a specific heat at constant pressure of
    cₚ = newScalar(ENTROPYperUnitMASS)
    set!(cₚ, 1.7e7)

    # Properties pertaining to the reference configuration κᵣ.

    # Fiber diameter is a random value drawn from a statisitcal distribution.
    # Statistics from: Sobin, Fung and Tremer, J. Appl. Phys., Vol. 64, 1988.
    # They found the square root of fiber diameter to distribute normally.
    μ = 0.952    # mean of the square root of fiber diameter in microns
    σ = 0.242    # its standard deviation in microns
    sqrtdia = rand(Normal(μ, σ))
    # bound the permissible variability
    while sqrtdia < μ-3σ
        sqrtdia = rand(Normal(μ, σ))
    end
    dia = sqrtdia^2
    # convert from a diameter in microns to a diameter in centimeters
    dia = dia / 10000.0
    # and finally, determine the fiber cross-sectional area in cm²
    area = pi * dia^2 / 4.0
    Aᵣ = newScalar(AREA)
    set!(Aᵣ, area)

    # thermodynamic causes
    θᵣ = newScalar(CENTIGRADE)  # nominal body temperature
    set!(θᵣ, 37.0)
    # Fiber reference stress is caused by the pleural pressure Pₚ.
    Pₚ = newScalar(BARYE)
    set!(Pₚ, 900.0)                   # nominally about 9 cm H₂O ≈ 900 barye
    Rₐ = 2.685 * Lᵣ / 2               # radius of alveolus
    V  = (4/3) * pi * Rₐ * Rₐ * Rₐ    # volume of alveolus as a sphere
    Vₐ = V * 2.785 / 4.189  # scaled: volume of dodecahedron / volume of sphere
    # Alveolar energy caused by nominal pleural pressure is
    Uₐ = Pₚ * Vₐ
    # There are 20 chords per alveolus, shared by three alveoli each, and
    # there are 12 membranes per alveolus, shared by two alveoli each, where
    # the energy per membrane is assumed to equal the energy per chord.
    # Hence, the internal energy per chord is
    Uc = Uₐ / (20/3 + 12/2)
    # 90% of which is taken to be carried by collagen, and 10% by elastin.
    σᵣ = 0.9 * Uc / (Lᵣ * Aᵣ)

    # thermodynamic effects
    ηᵣ = newScalar(ENTROPYperUnitMASS)
    set!(ηᵣ, 3.7e7)
    ϵᵣ = newScalar(DIMENSIONLESS)  # by definition, ϵᵣ = 0 in κᵣ

    # Properties pertaining to the initial configuration κ₀.

    # dimensional property
    A₀ = Aᵣ * Lᵣ / L₀  # fiber is isochoric

    # thermodynamic fields
    ϵ₀ = newScalar(DIMENSIONLESS)
    set!(ϵ₀, log(L₀/Lᵣ))
    Eᵣ = newScalar(MODULUS)  # modulus in reference configuration
    set!(Eᵣ, 5.0e7)          # the reference modulus for collagen
    σ₀ = σᵣ + Eᵣ * ϵ₀
    # Assume an isothermal elastic response from κᵣ → κ₀.
    θ₀ = deepcopy(θᵣ)
    η₀ = ηᵣ + α * (σ₀ - σᵣ) / ρ

    # Create a collagen fiber using these material parameters.
    return BioFiber(ruptured, ϵf, ρ, α, cₚ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end # _collagen

function _elastin(Lᵣ, L₀::PhyScalar)::BioFiber
    # The fiber is originally intact.
    ruptured = MBool(false)

    # Physical properties include a rupture strain of
    μ = 3.5    # mean value for the rupture strain
    σ = 0.5   # standard deviation of the rupture strain
    ϵr = rand(Normal(μ, σ))
    while ϵr < 1.0
        ϵr = rand(Normal(μ, σ))
    end
    ϵf = newScalar(DIMENSIONLESS)
    set!(ϵf, ϵr)
    # and the mass density for hydrated elastin in g/cm³ of
    ρ = newScalar(MASSDENSITY)
    set!(ρ, 1.31)

    # Thermal properties include a thermal strain in 1/°C of
    α = newScalar(THERMALSTRAIN)
    set!(α, 0.003)
    # and a specific heat at constant pressure of
    cₚ = newScalar(ENTROPYperUnitMASS)
    set!(cₚ, 4.2e7)

    # Properties pertaining to the reference configuration κᵣ.

    # Fiber diameter is a random value drawn from a statisitcal distribution.
    # Statistics from: Sobin, Fung and Tremer, J. Appl. Phys., Vol. 64, 1988.
    # They found the square root of fiber diameter to distribute normally.
    μ = 0.957    # mean of the square root of fiber diameter in microns
    σ = 0.239    # its standard deviation in microns
    sqrtdia = rand(Normal(μ, σ))
    # bound the permissible variability
    while sqrtdia < μ-3σ
        sqrtdia = rand(Normal(μ, σ))
    end
    dia = sqrtdia^2
    # convert from a diameter in microns to a diameter in centimeters
    dia = dia / 10000.0
    # and finally, determine the fiber cross-sectional area in cm²
    area = pi * dia^2 / 4.0
    Aᵣ = newScalar(AREA)
    set!(Aᵣ, area)

    # thermodynamic causes
    θᵣ = newScalar(CENTIGRADE)  # nominal body temperature
    set!(θᵣ, 37.0)
    # Fiber reference stress is caused by the pleural pressure Pₚ.
    Pₚ = newScalar(BARYE)
    set!(Pₚ, 900.0)                   # nominally about 9 cm H₂O ≈ 900 barye
    Rₐ = 2.685 * Lᵣ / 2               # radius of alveolus
    V  = (4/3) * pi * Rₐ * Rₐ * Rₐ    # volume of alveolus as a sphere
    Vₐ = V * 2.785 / 4.189  # scaled: volume of dodecahedron / volume of sphere
    # Alveolar energy caused by nominal pleural pressure is
    Uₐ = Pₚ * Vₐ
    # There are 20 chords per alveolus, shared by three alveoli each, and
    # there are 12 membranes per alveolus, shared by two alveoli each, where
    # the energy per membrane is assumed to equal the energy per chord.
    # Hence, the internal energy per chord is
    Uc = Uₐ / (20/3 + 12/2)
    # 90% of which is taken to be carried by collagen, and 10% by elastin.
    σᵣ = 0.1 * Uc / (Lᵣ * Aᵣ)

    # thermodynamic effects
    ηᵣ = newScalar(ENTROPYperUnitMASS)
    set!(ηᵣ, 3.4e7)
    ϵᵣ = newScalar(DIMENSIONLESS)  # by definition ϵᵣ = 0 in κᵣ

    # Properties pertaining to the initial configuration κ₀.

    # dimensional property
    A₀ = Aᵣ * Lᵣ / L₀  # fiber is isochoric

    # thermodynamic fields
    ϵ₀ = newScalar(DIMENSIONLESS)
    set!(ϵ₀, log(L₀/Lᵣ))
    Eᵣ = newScalar(MODULUS)  # modulus in reference configuration
    set!(Eᵣ, 1.6e6)          # the rubbery modulus for elastin
    σ₀ = σᵣ + Eᵣ * ϵ₀
    # Assume an isothermal elastic response from κᵣ → κ₀.
    θ₀ = deepcopy(θᵣ)
    η₀ = ηᵣ + α * (σ₀ - σᵣ) / ρ

    # Create an elastin fiber using these material parameters.
    return BioFiber(ruptured, ϵf, ρ, α, cₚ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end # _elastin

"""
    newAlveolarChord(Lᵣ, L₀::PhyScalar)::AlveolarChord

Creates a new instance of type AlveolarChord.  The first argument, Lᵣ, is the length of the chord in its reference configuration κᵣ, and the second argument, L₀ (L₀ ≥ Lᵣ), is its length in the initial configuration κ₀.  To the extent possible, parametric values are assigned values via distributions.
"""
function newAlveolarChord(Lᵣ, L₀::PhyScalar)::AlveolarChord
    # Verify the input values.
    if Lᵣ.u ≠ CENTIMETER
        msg = string("Argument Lᵣ must have units of CENTIMETER.")
        throw(ErrorException(msg))
    end
    if Sget(Lᵣ) < eps(Float32)
        msg = string("Argument Lᵣ must be positive valued.")
        throw(ErrorException(msg))
    end
    if L₀.u ≠ CENTIMETER
        msg = string("Argument L₀ must have units of CENTIMETER.")
        throw(ErrorException(msg))
    end
    if L₀ < Lᵣ
        L₀ = Lᵣ
    end
    # Create the fibers comprising the chord.
    fᶜ = _collagen(Lᵣ, L₀)
    fᵉ = _elastin(Lᵣ, L₀)

    # structural property of the chord:
    ruptured = MBool(false)

    # physical properties of the chord:
    ϕᶜ = fᶜ.Aᵣ / (fᶜ.Aᵣ + fᵉ.Aᵣ)
    ρ  = ϕᶜ * fᶜ.ρ + (1.0-ϕᶜ) * fᵉ.ρ

    # Chord properties in a reference configuration κᵣ:

    # dimensional property of the chord:
    Aᵣ = fᶜ.Aᵣ + fᵉ.Aᵣ

    # thermodynamic conjugate fields: causes
    θᵣ = ϕᶜ * fᶜ.θᵣ + (1.0-ϕᶜ) * fᵉ.θᵣ
    σᵣ = ϕᶜ * fᶜ.σᵣ + (1.0-ϕᶜ) * fᵉ.σᵣ

    # thermodynamic conjugate fields: effects
    ηᵣ = fᶜ.ηᵣ + fᵉ.ηᵣ
    ϵᵣ = newScalar(DIMENSIONLESS)  # by definition ϵᵣ = 0 in κᵣ

    # Fiber properties in the initial configurations κ₀:

    # dimensional properties:
    A₀ = fᶜ.A₀ + fᵉ.A₀

    # thermodynamic conjugate fields: causes
    θ₀ = ϕᶜ * fᶜ.θ₀ + (1.0-ϕᶜ) * fᵉ.θ₀
    σ₀ = ϕᶜ * fᶜ.σ₀ + (1.0-ϕᶜ) * fᵉ.σ₀

    # thermodynamic conjugate fields: effects
    η₀ = fᶜ.η₀ + fᵉ.η₀
    ϵ₀ = newScalar(DIMENSIONLESS)
    set!(ϵ₀, log(L₀/Lᵣ))

    return AlveolarChord(fᶜ, fᵉ, ruptured, ρ, ϕᶜ, Lᵣ, Aᵣ, θᵣ, σᵣ, ηᵣ, ϵᵣ, L₀, A₀, θ₀, σ₀, η₀, ϵ₀)
end # newAlveolarChord

end # module BioFibers
