function calculateObstacleXY!(obstacle, trackCoeff, xy_track::Array{Float64,2}, i::Int64, j::Int64)
	# all thes values are currently just used for plotting the semi axes of the obstacle
	#obstacle.index = obstacle.s_obstacle/trackCoeff.ds+1 
	s_over = ceil(Int64,obstacle.s_obstacle[i,j]/trackCoeff.ds)*trackCoeff.ds
	s_over_ind = convert(Int64,s_over/trackCoeff.ds+1)
	s_under = floor(Int64,obstacle.s_obstacle[i,j]/trackCoeff.ds)*trackCoeff.ds
	s_under_ind = convert(Int64,s_under/trackCoeff.ds+1)
	#if the obstacle is currently positioned  exactly at on eof the point on the grid we choose the next point for the secant to achieve approxiamted results
	if s_under_ind == s_over_ind
		s_over_ind = s_under_ind+1
		s_over = s_over+trackCoeff.ds
	end    
	if  s_over_ind > size(xy_track)[2]
			s_over_ind  = size(xy_track)[2]
			s_under_ind  = s_over_ind-1
	end
	weight1 = (obstacle.s_obstacle[i,j]-s_under)/trackCoeff.ds
	weight2 = (s_over-obstacle.s_obstacle[i,j])/trackCoeff.ds
	xy_coord_infront = [xy_track[1,s_over_ind] xy_track[2,s_over_ind]]
	xy_coord_before = [xy_track[1,s_under_ind] xy_track[2,s_under_ind]]
	#we take the weighted sums #in the case the object is at a node we just take the value from s_under which is the currnt postion
	obstacle.xy_vector[i,:,j] = xy_coord_infront *weight1+xy_coord_before*weight2
	x_secant = xy_track[1,s_over_ind] - xy_track[1,s_under_ind]   
	y_secant = xy_track[2,s_over_ind] - xy_track[2,s_under_ind]   
	secant_vec = [ x_secant y_secant]/norm([ x_secant y_secant])
	normal_vec = [-y_secant x_secant]/norm([-y_secant x_secant])
	vector_sy = vec(normal_vec*obstacle.sy_obstacle[i,j])# explicit conversion to vector is necesary because the dimensions in next equations would not coincide Array{FLoat64,1}vs Array{FLoat64,2}
	obstacle.xy_vector[i,:,j] = obstacle.xy_vector[i,:,j] + vector_sy

	#points to plot the ry semi axis of the ellipsis
	vector_ry = vec(normal_vec*obstacle.ry)
	obstacle.axis_y_up[i,:,j] = obstacle.xy_vector[i,:,j] + vector_ry
	obstacle.axis_y_down[i,:,j] = obstacle.xy_vector[i,:,j] - vector_ry

	#points to plot the rs semi axis of the ellipsis
	vector_rs = vec(secant_vec*obstacle.rs)
	obstacle.axis_s_up[i,:,j] = obstacle.xy_vector[i,:,j] + vector_rs
	obstacle.axis_s_down[i,:,j] = obstacle.xy_vector[i,:,j] - vector_rs

	nothing
end



function calculatePredictedXY(z_log::Array{Float64}, mpcParams, trackCoeff, xy_track::Array{Float64,2}, i::Int64, round::Int64)
	# all thes values are currently just used for plotting the predicted states
	s_over = zeros(mpcParams.N+1)
	s_under = zeros(mpcParams.N+1)
	s_over_ind = zeros(mpcParams.N+1)
	s_under_ind = zeros(mpcParams.N+1)

	xy_coord_infront = zeros(mpcParams.N+1,2)
	xy_coord_before = zeros(mpcParams.N+1,2)
	xy_pred = zeros(mpcParams.N+1,2)
	x_secant = zeros(mpcParams.N+1)
	y_secant = zeros(mpcParams.N+1)
	secant_vec = zeros(mpcParams.N+1,2)
	normal_vec = zeros(mpcParams.N+1,2)
	vector_sy = zeros(mpcParams.N+1,2)
	weight1 = zeros(mpcParams.N+1)
	weight2 = zeros(mpcParams.N+1)

	s_over = ceil(Int64,z_log[:,6,i]/trackCoeff.ds)*trackCoeff.ds
	s_under = floor(Int64,z_log[:,6,i]/trackCoeff.ds)*trackCoeff.ds
	
	for j =1:mpcParams.N+1
		if s_under[j] <0
			s_under[j] = 0
			#warn("predicted state with negative s \n in predicted step $j at iteration: $i of round: $round \n use s = 0 for plot \n s_pred = $(z_log[j,1,i])")
		end
		if s_over[j] <0
			s_over[j] = 0
			warn("predicted negative s \"big\": = $(z_log[j,6,i])")
		end

		#if over index zu gross dannn an anfang plotte
		s_over_ind[j] = convert(Int64,s_over[j]/trackCoeff.ds+1)
		s_under_ind[j] = convert(Int64,s_under[j]/trackCoeff.ds+1)
		if s_under_ind[j] == s_over_ind[j]
			s_over_ind[j] = s_under_ind[j]+1
			s_over[j] = s_over[j]+trackCoeff.ds
		end   
		if s_over_ind[j] > size(xy_track)[2]
			s_over_ind[j] = s_over_ind[j]%size(xy_track)[2]
		end
		if s_under_ind[j] > size(xy_track)[2]
			s_under_ind[j] = s_under_ind[j]%size(xy_track)[2]
		end
		s_over_ind  = convert(Array{Int64},s_over_ind)
		s_under_ind  = convert(Array{Int64},s_under_ind)
		weight1[j] = (z_log[j,6,i]-s_under[j])/trackCoeff.ds
		weight2[j] = (s_over[j]-z_log[j,6,i])/trackCoeff.ds
		#!!index line77
		xy_coord_infront[j,:] = [xy_track[1,s_over_ind[j]] xy_track[2,s_over_ind[j]]]
		xy_coord_before[j,:] = [xy_track[1,s_under_ind[j]] xy_track[2,s_under_ind[j]]]
		#we take the weighted sums #in the case the object is at a node we just take the value from s_under which is the currnt postion
		xy_pred[j,:] = xy_coord_infront[j,:] *weight1[j]+xy_coord_before[j,:]*weight2[j]
		x_secant[j] = xy_track[1,s_over_ind[j]] - xy_track[1,s_under_ind[j]]   
		y_secant[j] = xy_track[2,s_over_ind[j]] - xy_track[2,s_under_ind[j]]   
		secant_vec[j,:] = [ x_secant[j] y_secant[j]]/norm([ x_secant[j] y_secant[j]])
		normal_vec[j,:]= [-y_secant[j] x_secant[j]]/norm([-y_secant[j] x_secant[j]])
		vector_sy[j,:] = normal_vec[j,:]*z_log[j,5,i]# explicit conversion to vector is necesary because the dimensions in next equations would not coincide Array{FLoat64,1}vs Array{FLoat64,2}
		xy_pred[j,:] = xy_pred[j,:] + vector_sy[j,:]
	end

	return xy_pred
end

