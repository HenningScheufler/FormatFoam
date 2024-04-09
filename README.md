## install clang format

sudo apt-get install clang-format

###

refactor the entire project with the specified files

```
./formatOF.sh clangFormat_Proposed
```

```
./formatOF.sh clangFormat_Mozzilla
```


clang-format -i --style=file:<selectedformat> src/finiteVolume/fvMesh/fvMesh.C


## compile



navigate to the root
```
source etc/bashrc
```

and compile OpenFOAM

```
./Allwmake
```

## Caveat

Some of the macro have line break after the definition of the function resulting in a compilation error after formatting:
```
define addToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
```

By applying the patch, this issue can be resolved:

```
define addToRunTimeSelectionTable(baseType,thisType,argNames)                                                   \
```

The tested formatter are able to compile OpenFOAM
