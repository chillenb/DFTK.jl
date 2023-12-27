@testitem "Adaptive depth Anderson" setup=[TestCases] begin
using DFTK
using LinearAlgebra
(; silicon, aluminium) = TestCases.all_testcases

function test_addiis(testcase; temperature=0, Ecut=10, kgrid=[3, 3, 3])
    model = model_LDA(testcase.lattice, testcase.atoms, testcase.positions; temperature)
    basis = PlaneWaveBasis(model; kgrid, Ecut)
    tol   = 1e-10

    solver = scf_anderson_solver(; errorfactor=Inf, maxcond=1e6, m=100)
    scfres_rdiis = self_consistent_field(basis; tol, mixing=SimpleMixing(), solver)

    for errorfactor in (10, 100, 1000, 1e4, 1e5)
        println()
        println("#-- $errorfactor")
        println()
        solver = scf_anderson_solver(; errorfactor, maxcond=Inf, m=100)
        scfres = self_consistent_field(basis; tol, mixing=SimpleMixing(), solver)
        @test norm(scfres.ρ - scfres_rdiis.ρ) * sqrt(basis.dvol) < 10tol
    end
end

@testset "Silicon, no temp" begin
    test_addiis(silicon; temperature=0)
end
@testset "Aluminium, temp" begin
    test_addiis(aluminium; temperature=0.03)
end
end
