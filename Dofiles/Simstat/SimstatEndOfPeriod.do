			mat temp=(0.0,0.0,0.0)
			
			global Y = "salarycurr000"
			
			quietly su $Y if endspell!=1  // average salary among all staff who continue on to next year
			if r(N)==1 {
				mat temp=(r(N),r(mean),0.0)
			}
			else if r(N)>1 {
				mat temp=(r(N),r(mean),r(sd))
			}
			
			if $j==1 { //append parameters for this cell to the others
				mat eop=temp
			}
			else{
				mat eop=eop\temp
			}

			mat temp=(0.0,0.0,0.0)
			
			global Y = "salarycurr000"
			
			quietly su $Y if endspell==1  // average salary among all staff who exit
			if r(N)==1 {
				mat temp=(r(N),r(mean),0.0)
			}
			else if r(N)>1 {
				mat temp=(r(N),r(mean),r(sd))
			}
			
			if $j==1 { //append parameters for this cell to the others
				mat out=temp
			}
			else{
				mat out=out\temp
			}

						
			quietly su $Y if newspell!=1 // average salary among all incumbent staff 
			*su $Y if newspell!=1 // average salary among all incumbent staff 
			if r(N)==0 {
				mat temp=(0.0,9999,9999)
				}
			if r(N)==1 {
				mat temp=(r(N),r(mean),0.0)
			}
			else if r(N)>1 {
				mat temp=(r(N),r(mean),r(sd))
			}
			
			if $j==1 { //append parameters for this cell to the others
				mat inter=temp
			}
			else{
				mat inter=inter\temp
			}
