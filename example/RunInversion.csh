#!/bin.csh

set julia = /home/diba/julia/julia-1.7.3/bin/julia
set prg = /home/diba/REG1DMT/src/main.jl

set filedata = syndata.dat
set filestg = setting.dat
set fileout = run01.out

$julia $prg << END
$filedata
$filestg
$fileout
END
