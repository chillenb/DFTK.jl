@testitem "Adaptive depth Anderson" setup=[TestCases] begin
using DFTK
(; silicon, aluminium) = TestCases.all_testcases

function test_addiis(testcase; temperature=0, Ecut=10, kgrid=[3, 3, 3])
    model = model_LDA(testcase.lattice, testcase.atoms, testcase.positions; temperature)
    basis = PlaneWaveBasis(model; kgrid, Ecut)

    solver = scf_anderson_solver(; errorfactor=Inf, maxcond=1e6)
    scfres_rdiis = self_consistent_field(basis; tol=1e-4, mixing=SimpleMixing(), solver)
    println()
    println()

    for errorfactor in (10, 100, 1e4, 1e6)
        solver = scf_anderson_solver(; errorfactor, maxcond=Inf)
        scfres = self_consistent_field(basis; tol=1e-4, mixing=SimpleMixing(), solver)
        @test scfres.energies.total ≈ scfres_rdiis.energies.total atol=tol
        println()
        println()
    end
end

@testset "Silicon, no temp" begin
    test_addiis(silicon; temperature=0)
end
@testset "Aluminium, temp" begin
    test_addiis(aluminium; temperature=0.03)
end
end
