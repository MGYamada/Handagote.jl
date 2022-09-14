# Handagote.jl

[![Build Status](https://travis-ci.com/MGYamada/Handagote.jl.svg?branch=main)](https://travis-ci.com/MGYamada/Handagote.jl)
[![Coverage](https://codecov.io/gh/MGYamada/Handagote.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MGYamada/Handagote.jl)
[![Coverage](https://coveralls.io/repos/github/MGYamada/Handagote.jl/badge.svg?branch=main)](https://coveralls.io/github/MGYamada/Handagote.jl?branch=main)

Experimental forward-mode AD library for tensor networks.

Hadagote = Handai + Zygote. This is because the library was created
when the author was in Handai (Osaka University) in Japan.

Compatible with OMEinsum.jl, Zygote.jl, etc.
The implementation is following HyperDualNumbers.jl.
https://github.com/JuliaDiff/HyperDualNumbers.jl

## Installation

Type `]add https://github.com/MGYamada/Handagote.jl.git` from REPL.

## TODO

* Dual Number operation
* supporting op(::Dual, ::HyperDual)
* supporting CUDA

## Einsum

In order to avoid type piracy, the functionality of OMEinsum.jl is duplicated.
Currently only supports the `@dualein_str` notation `dualein"ij, jk, klm -> ilm"(A, B, C)`,
and `@hyperdualein_str` notation `hyperdualein"ij, jk, klm -> ilm"(A, B, C)`.

## Bug

`similar(A, Dual{Float64}, (2, 3))` does not work. Use `similar(A, Dual{Float64}, 2, 3)`, instead.

## Future

Zygote over Handagote implementation of MPRG is a temporal ad-hoc solution.
In the future, we will change everything to Diffractor.jl backend.

## License

MIT

## Author

* Masahiko G. Yamada