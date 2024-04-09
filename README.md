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

After formatting the runtime selectionTable OpenFOAM no longer compiles

clang-format -i --style=file:clangFormat_Proposed src/OpenFOAM/db/runTimeSelection/construction/runTimeSelectionTables.H