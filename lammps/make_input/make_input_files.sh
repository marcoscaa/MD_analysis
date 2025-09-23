####################
#User defined inputs
pos=$1       #User needs to input the filename of lammpstrj file
atype1=2     #Atom type used for angle distribution
atype2=2     #Atom type used for angle distribution
ntype1=2     # Number of atoms of type 1
ntype2=2     # Number of atoms of type 2
atypemsd=2   #Atom type used to compute MSD
nequil=0     #Number of initial steps not included for analysis
rcut="3.3"   #Radius cutoff to define 1st nearest neighbors for angle distribution
stride=1     #Run analysis at each $stride frames
nhistHb=1000 #Number of points in the histogram for H-bond distribution
nhistang=1000 #Number of points in the histogram for angle distribution
nhistrdf=5000 #Number of points in the histogram for radial distribution
nhistlsi=400 #Number of points in the histogram for LSI distribution
nhistzdf=100 #Numver of points in the histogram for number density distribution
dt="1.0"     #Timestep (in fs)  
rmax=40      #Maximum radius to consider for distribution of orientation 
center="0.0 0.0 0.0" #For OH orientation distribution
tcorr=1000   #Maximum length of time correlation functions
zoffset=0.0  #Offset along surface normal direction for ZDF
nlayers=1    #Number of layers (bins) along the Z direction
layers="0 1" #Specify the boundaries of layers along z (angstrom)

###########################################
#Variables obtained from the lammpstrj file
nat=`head $pos | awk '(NR==4){print $1}'`
nframes=`grep TIMESTEP $pos | wc -l | awk '{print $1}'`
nattype=`tail -n $nat $pos | awk -v"a=$atypemsd" 'BEGIN{c=0} ($2==a){c+=1} END{print c}'`
ntype=`tail -n $nat $pos | awk '{print $2}' | sort | tail -n 1`
a=`tail -n $((nat+9)) $pos | head -n 10 | awk '(NR==6){print $2-$1}'`
b=`tail -n $((nat+9)) $pos | head -n 10 | awk '(NR==7){print $2-$1}'`
c=`tail -n $((nat+9)) $pos | head -n 10 | awk '(NR==8){print $2-$1}'`

#######################
#Making the input files
cat << EOF > index_angle
$nat $nframes $nequil $nhistang
$atype1 $atype2 $rcut
EOF

cat << EOF > index_hb
$nat $nframes $nequil $stride $nhistHb $rmax 
$center
EOF

cat << EOF > index_msd
$nat $nattype $nframes $nequil $dt $stride $atypemsd $tcorr
EOF

cat << EOF > index_rdf
$nat $nframes $nequil $stride $nhistrdf
$a $b $c
$center
EOF

cat << EOF > index_orient
$nat $nframes $nequil $stride $nhistang $rmax
$center
EOF

cat << EOF > index_order
$nat $nframes $nequil $stride $nhistlsi
EOF

cat << EOF > index_zdf
$nat $ntype $nframes $nequil $stride $nhistzdf 3
$zoffset
EOF

cat << EOF > index_bond_angle
$nat $nframes $nequil $stride $nhistang
2 3
EOF

cat << EOF > index_vac
$nat $nframes $nequil $tcorr $dt $nlayers
$layers
EOF

cat << EOF > index_rcorr
$nat $ntype1 $ntype2 $nframes $nequil $dt $stride $nlayers
$atype1 $atype2 $tcorr
$layers
EOF

