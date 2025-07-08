companion_form_phi = function(Phi,d,p){

  # Main body
  comp_Phi <- array(0,dim=c((p*d),(p*d)))
  if (p == 1){
    comp_Phi = Phi[,,1]
  }
  else{
    for (i in 1:p){
      for (j in 1:p){
        if(i == 1){
          comp_Phi[1:d,((j-1)*d+1):(j*d)] <- Phi[,,j]
        }
        else if ( (i-1) == j & i != 1 ){
          comp_Phi[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- diag(1,d)
        }
      }
    }
  }

  return(comp_Phi)
}
