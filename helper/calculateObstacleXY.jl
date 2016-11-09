function calculateObstacleXY!(obstacle::Obstacle, trackCoeff::TrackCoeff, xy_track::Array{Float64,2}, i::Int64)
	# all thes values are currently just used for plotting the semi axes of the obstacle
	#obstacle.index = obstacle.s_obstacle/trackCoeff.ds+1 
	s_over = ceil(Int64,obstacle.s_obstacle[i]/trackCoeff.ds)*trackCoeff.ds
	s_over_ind = convert(Int64,s_over/trackCoeff.ds+1)
	s_under = floor(Int64,obstacle.s_obstacle[i]/trackCoeff.ds)*trackCoeff.ds
	s_under_ind = convert(Int64,s_under/trackCoeff.ds+1)
	#if the obstacle is currently positioned  exactly at on eof the point on the grid we choose the next point for the secant to achieve approxiamted results
	if s_under_ind == s_over_ind
		s_over_ind = s_under_ind+1
		s_over = s_over+trackCoeff.ds
	end           
	weight1 = (obstacle.s_obstacle[i]-s_under)/trackCoeff.ds
	weight2 = (s_over-obstacle.s_obstacle[i])/trackCoeff.ds
	xy_coord_infront = [xy_track[1,s_over_ind] xy_track[2,s_over_ind]]
	xy_coord_before = [xy_track[1,s_under_ind] xy_track[2,s_under_ind]]
	#we take the weighted sums #in the case the object is at a node we just take the value from s_under which is the currnt postion
	obstacle.xy_vector[i,:] = xy_coord_infront *weight1+xy_coord_before*weight2
	x_secant = xy_track[1,s_over_ind] - xy_track[1,s_under_ind]   
	y_secant = xy_track[2,s_over_ind] - xy_track[2,s_under_ind]   
	secant_vec = [ x_secant y_secant]/norm([ x_secant y_secant])
	normal_vec = [-y_secant x_secant]/norm([-y_secant x_secant])
	vector_sy = vec(normal_vec*obstacle.sy_obstacle[i])# explicit conversion to vector is necesary because the dimensions in next equations would not coincide Array{FLoat64,1}vs Array{FLoat64,2}
	obstacle.xy_vector[i,:] = obstacle.xy_vector[i,:] + vector_sy

	#points to plot the ry semi axis of the ellipsis
	vector_ry = vec(normal_vec*obstacle.ry)
	obstacle.axis_y_up[i,:] = obstacle.xy_vector[i,:] + vector_ry
	obstacle.axis_y_down[i,:] = obstacle.xy_vector[i,:] - vector_ry

	#points to plot the rs semi axis of the ellipsis
	vector_rs = vec(secant_vec*obstacle.rs)
	obstacle.axis_s_up[i,:] = obstacle.xy_vector[i,:] + vector_rs
	obstacle.axis_s_down[i,:] = obstacle.xy_vector[i,:] - vector_rs

	nothing
end

