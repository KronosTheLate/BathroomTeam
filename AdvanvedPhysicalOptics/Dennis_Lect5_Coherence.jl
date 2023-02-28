using Unitful
c0 = Unitful.c0
λ = (i=589.593u"nm", ii=588.996u"nm")
k(λ) = 2π/λ

T₊_spatial = 2π/(k(λ.i) + k(λ.ii)) |> u"μm"
T₋_spatial = 2π/(-k(λ.i) + k(λ.ii)) |> u"mm"

N_fast_oscillations_per_slow_oscillation = T₋_spatial/T₊_spatial |> NoUnits

T₊_temporal = T₊_spatial/c0 |> u"fs"
T₋_temporal = T₋_spatial/c0 |> u"fs"