# Handagote.jl

Compatible with OMEinsum.jl, Zygote.jl, etc.

## TODO

hack broadcast / supporting CUDA

## Einsum

In order to avoid type piracy, the functionality of OMEinsum.jl is duplicated.
Currently only supports the `@dualein_str` notation `dualein"ij, jk, klm -> ilm"(A, B, C)`,
and `@hyperdualein_str` notation `hyperdualein"ij, jk, klm -> ilm"(A, B, C)`