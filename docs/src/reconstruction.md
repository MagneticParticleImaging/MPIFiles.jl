# Reconstruction Results

MDF files can also contain reconstruction results instead of measurement data.
The low level results can be retrieved using the [Low Level Interface](@ref), e.g.
[`recoData`](@ref), [`recoFov`](@ref), [`recoFovCenter`](@ref), [`recoSize`](@ref),
[`recoOrder`](@ref) and [`recoPositions`](@ref).

Instead, one can also combine these data into an `ImageMetadata` object from the
[ImageMetadata.jl](https://github.com/JuliaImages/ImageMetadata.jl) package.
The `ImageMeta` object also pulls all relevant metadata from an MDF
such that the image can be stored again. These two functions are especially relevant when using the package
[MPIReco.jl](https://github.com/MagneticParticleImaging/MPIReco.jl).

```@docs
loadRecoData
saveRecoData
```
