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