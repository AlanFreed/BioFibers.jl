module testBioFibers

import
    PhysicalSystemsOfUnits:
        CGS,
        CENTIMETER

import
    PhysicalScalars:
        PhysicalScalar,
        newPhysicalScalar,
        set!

include("../src/BioFibers.jl")

export
    run

function run()
    Lᵣ = newPhysicalScalar(CENTIMETER)
    set!(Lᵣ, 65.0e-4)
    L₀ = newPhysicalScalar(CENTIMETER)
    set!(L₀, 66.0e-4)
    ac = BioFibers.newAlveolarChord(Lᵣ, L₀)
    print(BioFibers.toString(ac))
    print("\nA shallow copy is:\n\n")
    print(BioFibers.toString(BioFibers.copy(ac)))
    print("\nA deep copy is:\n\n")
    print(BioFibers.toString(BioFibers.deepcopy(ac)))
end  # run

end  # testBioFibers