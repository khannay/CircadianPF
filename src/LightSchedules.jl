#=
Collection of simple light schedules
=#



module LightSchedules

function RegularLight(t; Intensity=150.0, PP=16.0, period=24.0)
	#=
	Define a basic regular light schedule, can change the light intensity, PP
	the photoperiod and the period of the light forcing 24.0
	=#
	val=0.5*tanh(100*sin(2π/period*t))+0.5
    return(Intensity*val)
end

function RegularLightSimple(t; Intensity=150.0, wakeUp=8.0, workday=16.0)
	#=
	Define a basic light schedule with Intensity=150.0, wakeUp=8.0 and a
	workday=16.0 by default
	=#

	s=rem(t,24.0)-wakeUp
	if s<0.0
		s+=24.0
	end
	val=0.5*tanh(100*s)-0.5*tanh(100*(s-workday))
	return(Intensity*val)
end

function ShiftWorkLight(t; dayson=5, daysoff=2)
	t=rem(t, (dayson+daysoff)*24.0) #repeats
	if t<=24*dayson
		return RegularLightSimple(t, wakeUp=16.0, workday=16.0)
	else
		return RegularLightSimple(t, wakeUp=9.0, workday=16.0)
	end
end #shift worker

function SocialJetLag(t; weekdayWake=7.0, weekdayBed=24.0, weekendWake=11.0, weekendBed=2.0)
	t=rem(t, 7*24.0)
	if t<=24*5
		duration=rem(weekdayBed-weekdayWake,24.0)
		if duration<0.0
			duration+=24.0
		end
		return RegularLightSimple(t, wakeUp=weekdayWake, workday=duration)
	else
		duration=rem(weekendBed-weekendWake, 24.0)
		if duration<0.0
			duration+=24.0
		end
		return RegularLightSimple(t, wakeUp=weekendWake, workday=duration)
	end

end #SocialJetLag

function SlamShift(t; shift=8.0, Intensity=150.0, beforeDays=10)
	t=rem(t, 40*24)

	if t<= 24*beforeDays
		return RegularLightSimple(t, Intensity=Intensity, wakeUp=8.0, workday=16.0)
	else
		newVal=rem(8.0+shift, 24.0)
		if newVal<0
			newVal+=24.0
		end
		return RegularLightSimple(t, Intensity=Intensity, wakeUp=newVal, workday=16.0)
	end

end


end #module LightSchedules
