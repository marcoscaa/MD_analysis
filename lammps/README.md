# Post-processing codes for lammps

To compile the codes, first 

```
cd src
make gnu
cd ..
```

and then you can go to any other folder and

```
make gnu
```

or just ```make```. Look at the Makefile to make sure you are using the GNU compiler

Each executable takes as input first the lammpstrj file and then a code-specific input file.
You can generate this code-specific input file using the bash script make_input/make_input_files.sh.
Look inside this bash script and change any variable you need.
