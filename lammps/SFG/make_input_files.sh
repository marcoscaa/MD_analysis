####################
#User defined inputs
pos=$1       #User needs to input the filename of lammpstrj file
atypeO=2     # Index type of O atoms
atypeH=3     # Index type of H atoms
nequil=0     #Number of initial steps not included for analysis
stride=1     #Run analysis at each $stride frames
dt="0.001"   #Timestep (in ps)  
tcorr=2000   #Maximum length of time correlation functions
zoffset=0.0  #Offset along surface normal direction for ZDF
nlayers=1    #Number of layers (bins) along the Z direction
layers="0 1" #Specify the boundaries of layers along z (angstrom)

###########################################
#Variables obtained from the lammpstrj file
nat=`head $pos | awk '(NR==4){print $1}'`
nframes=`grep TIMESTEP $pos | wc -l | awk '{print $1}'`
nwater=`tail -n $nat $pos | awk -v "a=$atypeH" 'BEGIN(c=0); ($2==a){c+=1}; END{print a/2}' `

#######################
#Making the input files
cat << EOF > input_sfg
$nat $nframes $nequil 
$nwater $nlayers $dt
$layers
EOF
