
# La fonction recuit retourne le minimum de funk rencontré au-cours de l'exploration.
# x pointe alors vers le vecteur réalisant ce minimum (x initial étant le point de départ pour l'algorithme de simplexe)
# T_init est la température de départ, T_end celle de fin, alpha est le coefficient de refroidissement T(n+1)=alpha*T(n)
# L est la longueur du palier (pendant combien de cycles on reste à la même température)
# */
# // char_length est un vecteur de longueurs caractéristiques (pour chacune des composantes

function recuit_nm(funk::Function, x::Vector, ftol::Float64, T_init::Float64, T_end::Float64, alpha::Float64, L::Int64)
	ndim=length(x)
	iter=0;
	temp=T_init;

	# Definition du SImplex initial
	char_length=ones(ndim)  
	mat=Array(Float64, ndim+1, ndim)
	for i=1:ndim+1
		for j=1:ndim
			mat[i,j]=x[j]
		end
		if(i>1)
			mat[i,i-1]+=char_length[i-1];
		end
	end

	pb=zeros(ndim);
	y=zeros(ndim+1)
	
	#ybest
	yb=1000000.;
	ibest=1
	for i=1:ndim+1
		y[i]=funk(vec(mat[i,:]))
		if(y[i]<yb) 
			yb=y[i]
			ibest=i
		end
	end
	pb=copy(vec(mat[ibest,:]))

	while (iter<=0 && temp>T_end)
		println("\n********T=",temp);
		if(iter<=0)
			iter=L
			yb,iter=amebsa!(mat, y, ndim, pb, yb, ftol, funk, iter, temp)
		end
		println("f(pb)=",funk(pb))
		temp*=alpha
	end

	return pb, yb

end


# Multidimensional minimization of the function funk(x) where x[1..ndim] is a vector in
# ndim dimensions, by simulated annealing combined with the downhill simplex method of Nelder
# and Mead. The input matrix p[0..ndim][0..ndim-1] has ndim+1 rows, each an ndimdimensional
# vector which is a vertex of the starting simplex. Also input are the following: the
# vector y[0..ndim], whose components must be pre-initialized to the values of funk evaluated
# at the ndim+1 vertices (rows) of p; ftol, the fractional convergence tolerance to be
# achieved in the function value for an early return; iter, and temptr. The routine makes iter
# function evaluations at an annealing temperature temptr, then returns. You should then de-
# crease temptr according to your annealing schedule, reset iter, and call the routine again
# (leaving other arguments unaltered between calls). If iter is returned with a positive value,
# then early convergence and return occurred. If you initialize yb to a very large value on the first
# call, then yb and pb[0..ndim-1] will subsequently return the best function value and point ever
# encountered (even if it is no longer a point in the simplex).


function amebsa!(p::Array{Float64}, y::Vector, ndim::Int64, pb::Vector, yb::Float64, ftol::Float64, funk::Function, iter::Int64, temptr::Float64)
	mpts=ndim+1
	psum=zeros(ndim)
	tt = -temptr;
	for n=1:ndim
		sum=0.0
		for m=1:mpts
			sum += p[m,n]
		end
		psum[n]=sum
	end
	itmax=iter
	for it=1:itmax
		ilo=1 #Determine which point is the highest (worst), next-highest, and lowest (best). 
		ihi=2
		ynhi=ylo=y[1]+tt*log(rand()) #Whenever we “look at” a vertex, it gets a random thermal .uctuation. 
		yhi=y[2]+tt*log(rand())
		if (ylo > yhi) 
			ihi=1
			ilo=2
			ynhi=yhi
			yhi=ylo
			ylo=ynhi
		end
		for i=3:mpts
			/#Loop over the points in the simplex.
			yt=y[i]+tt*log(rand()) #More thermal .uctuations.
			if (yt <= ylo) 
				ilo=i;
				ylo=yt;
			end
			if (yt > yhi) 
				ynhi=yhi;
				ihi=i;
				yhi=yt;
			elseif (yt > ynhi) 
				ynhi=yt;
			end
		end
		rtol=2.0*abs(yhi-ylo)/(abs(yhi)+abs(ylo));
		#Compute the fractional range from highest to lowest and return if satisfactory.
		if (rtol < ftol || iter < 0) 
			#If returning, put best point and value in slot 1. 
			swap=y[1];
			y[1]=y[ilo];
			y[ilo]=swap;
			for n=1:ndim 
				swap=p[1,n];
				p[1,n]=p[ilo,n];
				p[ilo,n]=swap;
			end
			break
		end
		iter -= 2;
		#Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex across from the high point, i.e., re.ect the simplex from the high point.
		ytry, yb, yhi=amotsa!(p,y,psum,ndim,pb,yb,funk,ihi,yhi,-1.0, tt)
		if (ytry <= ylo) 
			#Gives a result better than the best point, so try an additional extrapolation by afactor of 2.#			println("yb , yhi avant:", yb, " , ", yhi)
			ytry, yb, yhi=amotsa!(p,y,psum,ndim,pb,yb,funk,ihi,yhi,2.0, tt);#			println("yb , yhi après:", yb, " , ", yhi)
		elseif (ytry >= ynhi) 
			#The reflected point is worse than the second-highest, so look for an intermediate 
			ysave=yhi;
			ytry, yb, yhi=amotsa!(p,y,psum,ndim,pb,yb,funk,ihi,yhi,0.5, tt);
			if (ytry >= ysave) 
			# Can’t seem to get rid of that high point. better contract around the lowest (best) point.
				for i=1:mpts
					if (i != ilo) 
						for j=1:ndim
							psum[j]=0.5*(p[i,j]+p[ilo,j]);
							p[i,j]=psum[j];
						end
						y[i]=funk(psum);
					end
				end
				iter -= ndim;
				#Recompute psum.
				for n=1:ndim
					sum=0.0
					for m=1:mpts
						sum += p[m,n]
					end
					psum[n]=sum
				end
			end
		else 
			iter+=1 #Correct the evaluation count.
		end
	end
	return yb, iter
end


#Extrapolates by a factor fac through the face of the simplex across from the high point, tries
#it, and replaces the high point if the new point is better.

function amotsa!(p::Array{Float64}, y::Vector, psum::Vector, ndim::Int64, pb::Vector, yb::Float64, funk::Function, ihi::Int64, yhi::Float64, fac::Float64, tt::Float64)
	ptry=zeros(ndim)
	fac1=(1.0-fac)/ndim
	fac2=fac1-fac
	for i=1:ndim
		ptry[i]=psum[i]*fac1-p[ihi,i]*fac2
	end
	ytry=funk(ptry)
	if (ytry <= yb) 
		#Save the best-ever.
		for j=1:ndim
			pb[j]=ptry[j]
		end
		yb=ytry
	end
	yflu=ytry-tt*log(rand()) 
	# We added a thermal fluctuation to all the current vertices, but we subtract it here, so as to give
	# the simplex a thermal Brownian motion: It likes to accept any suggested change.
	if (yflu < yhi) 
		y[ihi]=ytry
		yhi=yflu
		for j=1:ndim 
			psum[j] += ptry[j]-p[ihi,j]
			p[ihi,j]=ptry[j]
		end
	end
	newyb=yb
	newyhi=yhi
	return yflu, newyb, newyhi
end






# Recuit simulé simple: on utilise un pas de recherche selon chacune des composantes -> discrétisation

function recuit0(func::Function, x::Vector, xinf::Vector, xsup::Vector, pas::Vector, T_init::Float64, T_end::Float64, alpha::Float64, LMIN::Int64)

	ndim=length(x)
	iter=0
	temp=T_init

	y=func(x)
	y_min=y
	y_min_newton=y

	nbpal=0 #nb de pallier sans nouveaux résultats
	newmin=false
	L=LMIN
	LMAX = 10 * LMIN
	NBPALMAX=10 #Nb de pallier max sans nouveau min avant de diminuer la longueur de pallier
#	teta=Vector(Float64, 5)
#	teta_voisin=Vector(Float64, ndim)
#	x_newton=Vector(Float64, ndim)
	teta=copy(x)
	teta_voisin=copy(x)
	 # for i=1:ndim
		# teta[i] = x[i]
	# end


	while (iter<=0 && temp>T_end)
		newmin=0
		println("********T=",temp)
		for j=1:L
			teta_voisin=copy(teta)
			k = floor(rand()*ndim)+1  #on choisit au hasard un des paramètres
			l = sign(2*rand() - 1)	 #l vaut -1 ou 1
			if (teta_voisin[k]-pas[k])<=xinf[k]			
				teta_voisin[k]=xinf[k]+10*pas[k]; #on suppose que l'intervalle est suffisamment grand pour que la valeur à trouver ne soit pas sur le bord
											#ATTENTION : pose pb s'il y a moins de 10 pas dans chaque intervalle !!
			elseif (teta_voisin[k]+pas[k])>=xsup[k]
				teta_voisin[k] = xsup[k]-10*pas[k]
	#			if ( teta_voisin[k]-pas[k]<inf[k]  ||  teta_voisin[k]+pas[k]>sup[k] )
	#				teta_voisin[k]= (inf[k]+ sup[k])/2 ; //on se replace au milieu de l'intervalle
													 #le pb est que du coup, on s'éloigne de la solution...A VOIR !
			else
				teta_voisin[k] += l*pas[k]
			end
			
			y_voisin=func(teta_voisin)

			if(y_voisin<y)
				teta=copy(teta_voisin)
				y=y_voisin
			else
				proba=rand()	#donne un nombre au hasard entre 0 et 1.
				metro_proba = exp(-(y_voisin-y)/temp)	
			# if(temp == T_init)
				# println("\n Initial proba : ",metro_proba);
			# end
				if(proba <= metro_proba)
					teta=copy(teta_voisin)
					y=y_voisin
				end
			end
			
			if(y<y_min)
				println("Nouvelle estimation, d'erreur : ", y)
				for i=1:ndim
					x[i]=teta[i]	  #on garde en mémoire le plus petit de tous, qqsoit le palier
					println(x[i],"\t")
				end
				#lancer un algo de newton 
				# res=optimize(f4, df4!, teta, method = :bfgs)
				# x=copy(res.minimum)
				# teta=copy(res.minimum)
				# println("\n")
				y_min=y
				newmin=1
				nbpal=0
				L=LMAX
			end
		end
		if(newmin!=0)
			nbpal+=1
			if(nbpal>=NBPALMAX)
				nbpal=0
				L=max(LMIN, L/2)
			end
		end

		temp*=alpha
	end
	return y_min
end
