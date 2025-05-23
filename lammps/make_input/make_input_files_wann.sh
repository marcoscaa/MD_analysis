####################
#User defined inputs
pos=$1      #User needs to input the filename of lammpstrj file
dipole=$2   #User needs to input the trajectory containing dipoles
typeO=2     #Atom type used for O 
typeH=3     #Atom type used for H
typeC=1     #Atom type used for C
nequi=1000  #Number of initial steps to be discarded
temp=300    #Temperature in K
dt="1.0"    #Time difference between consecutive trajectory frames

###########################################
#Variables obtained from the lammpstrj file
nat=`head $pos | awk '(NR==4){print $1}'`
nwann=`head $dipole | awk '(NR==4){print $1}'`
nframes=`grep TIMESTEP $pos | wc -l | awk '{print $1}'`

#######################
#Making the input files
cat << EOF > index_wann
$nat $nwann $nframes $temp $dt 
$typeO $typeH $typeC
EOF
