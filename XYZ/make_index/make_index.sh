function make_index(){
  traj=$1
  #Make index file
  nat=`head -n 1 $traj | awk '{print $1}'`
  nframes=`grep Lattice $traj | wc -l | awk '{print $1}'`
  nattype=`head -n $((nat+2)) $traj | awk 'BEGIN{c=0} ($1=="N"){c+=1} END{print c}'`
  nequil=$((nframes/10))
  stride=1
  nbins=300
  rcut="6.0"

  cat << EOF > index_zdf
$nat $nframes $nequil $stride $nbins 3 
0.0 
EOF

  echo "${nat} ${nframes} ${nequil}" > index_oh
  #echo "${nat} ${nattype} ${nframes} 0 0.10 1 10000" > index_msd
  #echo "${nat} ${nframes} ${nequil} ${stride}" > index_hb

}

make_index $1
