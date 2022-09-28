# The Hessian of P -> E(P) (E being the energy) is Ω+K, where Ω and K are
# defined below (cf. [1] for more details).
#
# In particular, by linearizing the Kohn-Sham equations, we have
#
# δP = -(Ω+K)⁻¹ δH
#
# which can be solved either directly (solve_ΩplusK) or by a splitting method
# (solve_ΩplusK_split), the latter being preferable as it is well defined for
# both metals and insulators. Solving this equation is necessary to compute
# response properties as well as AD chain rules (Ω+K being self-adjoint, it can
# be used for both forward and reverse mode).
#
# [1] Eric Cancès, Gaspard Kemlin, Antoine Levitt. Convergence analysis of
#     direct minimization and self-consistent iterations.
#     SIAM Journal on Matrix Analysis and Applications
#     https://doi.org/10.1137/20M1332864
#
# TODO find better names for solve_ΩplusK and solve_ΩplusK_split
#

"""
    apply_Ω(δψ, ψ, H::Hamiltonian, Λ)

Compute the application of Ω defined at ψ to δψ. H is the Hamiltonian computed
from ψ and Λ is the set of Rayleigh coefficients ψk' * Hk * ψk at each k-point.
"""
function apply_Ω(δψ, ψ, H::Hamiltonian, Λ)
    δψ = proj_tangent(δψ, ψ)
    Ωδψ = [H.blocks[ik] * δψk - δψk * Λ[ik] for (ik, δψk) in enumerate(δψ)]
    proj_tangent!(Ωδψ, ψ)
end

"""
    apply_K(basis::PlaneWaveBasis, δψ, ψ, ρ, occupation)

Compute the application of K defined at ψ to δψ. ρ is the density issued from ψ.
δψ also generates a δρ, computed with `compute_δρ`.
"""
@views function apply_K(basis::PlaneWaveBasis, δψ, ψ, ρ, occupation)
    δψ = proj_tangent(δψ, ψ)
    δρ = compute_δρ(basis, ψ, δψ, occupation)
    δV = apply_kernel(basis, δρ; ρ=ρ)

    Kδψ = map(enumerate(ψ)) do (ik, ψk)
        kpt = basis.kpoints[ik]
        δVψk = similar(ψk)

        for n = 1:size(ψk, 2)
            ψnk_real = ifft(basis, kpt, ψk[:, n])
            δVψnk_real = δV[:, :, :, kpt.spin] .* ψnk_real
            δVψk[:, n] = fft(basis, kpt, δVψnk_real)
        end
        δVψk
    end
    # ensure projection onto the tangent space
    proj_tangent!(Kδψ, ψ)
end

"""
    solve_ΩplusK(basis::PlaneWaveBasis{T}, ψ, res, occupation;
                 tol=1e-10, verbose=false) where {T}

Return δψ where (Ω+K) δψ = rhs
"""
function solve_ΩplusK(basis::PlaneWaveBasis{T}, ψ, rhs, occupation;
                      tol=1e-10, verbose=false) where {T}
    @assert mpi_nprocs() == 1  # Distributed implementation not yet available
    filled_occ = filled_occupation(basis.model)
    # for now, all orbitals have to be fully occupied -> need to strip them beforehand
    @assert all(all(occ_k .== filled_occ) for occ_k in occupation)

    # compute quantites at the point which define the tangent space
    ρ = compute_density(basis, ψ, occupation)
    _, H = energy_hamiltonian(basis, ψ, occupation; ρ=ρ)

    pack(ψ) = reinterpret_real(pack_ψ(ψ))
    unpack(x) = unpack_ψ(reinterpret_complex(x), size.(ψ))
    unsafe_unpack(x) = unsafe_unpack_ψ(reinterpret_complex(x), size.(ψ))

    # project rhs on the tangent space before starting
    proj_tangent!(rhs, ψ)
    rhs_pack = pack(rhs)

    # preconditioner
    Pks = [PreconditionerTPA(basis, kpt) for kpt in basis.kpoints]
    for ik = 1:length(Pks)
        precondprep!(Pks[ik], ψ[ik])
    end
    function f_ldiv!(x, y)
        δψ = unpack(y)
        proj_tangent!(δψ, ψ)
        Pδψ = [ Pks[ik] \ δψk for (ik, δψk) in enumerate(δψ)]
        proj_tangent!(Pδψ, ψ)
        x .= pack(Pδψ)
    end

    # Rayleigh-coefficients
    Λ = [ψk'Hψk for (ψk, Hψk) in zip(ψ, H * ψ)]

    # mapping of the linear system on the tangent space
    function ΩpK(x)
        δψ = unsafe_unpack(x)
        Kδψ = apply_K(basis, δψ, ψ, ρ, occupation)
        Ωδψ = apply_Ω(δψ, ψ, H, Λ)
        pack(Ωδψ + Kδψ)
    end
    J = LinearMap{T}(ΩpK, size(rhs_pack, 1))

    # solve (Ω+K) δψ = rhs on the tangent space with CG
    δψ, history = cg(J, rhs_pack, Pl=FunctionPreconditioner(f_ldiv!),
                  reltol=0, abstol=tol, verbose=verbose, log=true)

    (; δψ=unpack(δψ), history)
end


"""
Solve the problem `(Ω+K) δψ = rhs` using a split algorithm, where `rhs` is typically
`-δHextψ` (the negative matvec of an external perturbation with the SCF orbitals `ψ`) and
`δψ` is the corresponding total variation in the orbitals `ψ`. Additionally returns:
    - `δρ`:  Total variation in density)
    - `δHψ`: Total variation in Hamiltonian applied to orbitals
    - `δeigenvalues`: Total variation in eigenvalues
    - `δVind`: Change in potential induced by `δρ` (the term needed on top of `δHextψ`
      to get `δHψ`).
"""
@timing function solve_ΩplusK_split(ham::Hamiltonian, ρ::AbstractArray{T}, ψ, occupation, εF,
                                    eigenvalues, rhs;
                                    tol=1e-8, verbose=false, occupation_threshold,
                                    tol_sternheimer_max=max(tol, 1e-2),  # Least accuracy
                                    λmin_epsilon=one(T),  # Estimated smallest eigenvalue of
                                                          # the dielectric operator
                                    kwargs...) where {T}
    # Using χ04P = -Ω^-1, E extension operator (2P->4P) and R restriction operator:
    # (Ω+K)^-1 =  Ω^-1 ( 1 -   K   (1 + Ω^-1 K  )^-1    Ω^-1  )
    #          = -χ04P ( 1 -   K   (1 - χ04P K  )^-1   (-χ04P))
    #          =  χ04P (-1 + E K2P (1 - χ02P K2P)^-1 R (-χ04P))
    # where χ02P = R χ04P E and K2P = R K E
    basis = ham.basis
    @assert size(rhs[1]) == size(ψ[1])  # Assume the same number of bands in ψ and rhs
    @assert tol_sternheimer_max ≥ tol

    if verbose
        println("n     Arnoldi   Residual    CG     Comment")
        println("---   -------   --------    -----  -------")
    end

    # compute δρ0 (ignoring interactions)
    χ0res = apply_χ0_4P(ham, ψ, occupation, εF, eigenvalues, -rhs;      # = -χ04P * rhs
                        reltol=0, abstol=tol, occupation_threshold, kwargs...)
    δρ0 = compute_δρ(basis, ψ, χ0res.δψ, occupation, χ0res.δoccupation)
    if verbose && mpi_master()
        avg_iter = mpi_mean(mean(sum, χ0res.iterations), basis.comm_kpts)
        @printf "%s%6.1f  Non-interacting response\n" " "^27 avg_iter
    end

    # Solve for total δρ (Dyson equation)
    # Adaptive tolerance following Simoncini and Szyld (2007), DOI 10.1002/nla.499
    χ0tol   = Ref(min(tol / 10, λmin_epsilon * tol / norm(δρ0)))
    χ0niter = Ref(0.0)  # Extract average number of iterations

    # Dielectric operator
    pack(δρ)   = vec(δρ)
    unpack(δρ) = reshape(δρ, size(ρ))
    epsilon = LinearMap{T}(prod(size(δρ0))) do δρ
        δρ = unpack(δρ)
        δV = apply_kernel(basis, δρ; ρ)
        χ0res = apply_χ0(ham, ψ, occupation, εF, eigenvalues, δV;
                        occupation_threshold, abstol=χ0tol[], reltol=0, kwargs...)

        # Addition because operator might be called multiple times between outputs
        χ0niter[] += mpi_mean(mean(sum, χ0res.iterations), basis.comm_kpts)

        pack(δρ - χ0res.δρ)
    end

    # Note: δρ is updated in-place by gmres_iterable!
    δρ = zero(δρ0)
    iter = IterativeSolvers.gmres_iterable!(pack(δρ), epsilon, pack(δρ0);
                                            reltol=0, abstol=tol, initially_zero=true)
    n_iter = 0
    for residual in iter
        n_iter += 1
        if verbose && mpi_master() && !IterativeSolvers.converged(iter)
            @printf "% 3d    % 6d   %8.2f   %6.1f  " n_iter iter.k log10(residual) χ0niter[]
            println(n_iter == 1 ? "Solving Dyson equation" : "")
        end
        χ0niter[] = 0.0

        will_restart = iter.k == iter.restart
        if will_restart  # Whenever we restart we need a more accurate solve
            χ0tol[] = tol / 10
        else
            χ0tol[] = min(tol_sternheimer_max, λmin_epsilon * tol / residual)
        end
    end

    # Compute total change in Hamiltonian applied to ψ
    δVind = apply_kernel(basis, δρ; ρ)  # Change in potential induced by δρ
    δHψ = @views map(basis.kpoints, ψ, rhs) do kpt, ψk, rhsk
        δVindψk = RealSpaceMultiplication(basis, kpt, δVind[:, :, :, kpt.spin]) * ψk
        δVindψk - rhsk
    end

    # Compute total change in eigenvalues
    δeigenvalues = map(ψ, δHψ) do ψk, δHψk
        map(eachcol(ψk), eachcol(δHψk)) do ψnk, δHψnk
            real(dot(ψnk, δHψnk))  # δε_{nk} = <ψnk | δH | ψnk>
        end
    end

    χ0res = apply_χ0_4P(ham, ψ, occupation, εF, eigenvalues, δHψ;
                        occupation_threshold, abstol=tol, reltol=0, kwargs...)
    if verbose && mpi_master()
        avg_iter = mpi_mean(mean(sum, χ0res.iterations), basis.comm_kpts)
        @printf "%s%6.1f  Interacting response\n" " "^27 avg_iter
    end

    (; χ0res.δψ, δρ, δHψ, δVind, δeigenvalues, χ0res.δoccupation, χ0res.δεF)
end

function solve_ΩplusK_split(basis::PlaneWaveBasis, ψ, rhs, occupation; kwargs...)
    ρ = compute_density(basis, ψ, occupation)
    _, H = energy_hamiltonian(basis, ψ, occupation; ρ)
    eigenvalues = [real.(eigvals(ψk'Hψk)) for (ψk, Hψk) in zip(ψ, H * ψ)]
    occupation, εF = compute_occupation(basis, eigenvalues)

    solve_ΩplusK_split(H, ρ, ψ, occupation, εF, eigenvalues, rhs; kwargs...)
end

function solve_ΩplusK_split(scfres::NamedTuple, rhs; kwargs...)
    solve_ΩplusK_split(scfres.ham, scfres.ρ, scfres.ψ, scfres.occupation,
                       scfres.εF, scfres.eigenvalues, rhs;
                       scfres.occupation_threshold, kwargs...)
end
