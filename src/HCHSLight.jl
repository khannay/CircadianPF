module HCHSLight
#=
This Module should readin and return a interpolated Light Schedule for HCHS data
=#

using CSV, DataFrames, Dierckx, DataFramesMeta, RollingFunctions


function readData(filename)
	# materialize a csv file as a DataFrame
	df = CSV.File(filename) |> DataFrame
	df.Lux=df.whitelight+df.bluelight+df.redlight+df.greenlight

	#Make a day counter for the system
	day=[0]
	count=0
	for i in 2:length(df.dayofweek)
		if df.dayofweek[i] != df.dayofweek[i-1]
			count+=1
		end
		push!(day, count)
	end

	df.day=day

	#Create an index to count the hours of the data (UNITS=Hours)
	TimeTotal=[]
	TimeCount=[]
	for i in 1:length(df.time)
		a=sum(map(x-> parse(Float64,x), split(df.time[i], ":")).*[1.0, 1/60, 1/3600.0])
		b=a+df.day[i]*24.0
		push!(TimeTotal, b)
		push!(TimeCount,a)
	end

	df.TimeTotal=TimeTotal
	df.TimeCount=TimeCount

	df2=@select(df, :TimeTotal, :Lux)
	df2=dropmissing(df2, disallowmissing=true)
	df2.Lux=runmean(df2.Lux, 5) #smooth the light data a bit using running mean

	#Need to add trim and periodic extend components
	xvals=df2.TimeTotal
	yvals=df2.Lux

	yvals=yvals[xvals .>= 24.0]
	xvals=xvals[xvals .>= 24.0]

	#Chop hangover from the end
	numDays=xvals[end]/24.0
	finalTime=floor(numDays)*24.0
	#Exclude the final day
	yvals=yvals[xvals .<= finalTime]
	xvals=xvals[xvals .<= finalTime]


	#start at zero
	xvals=xvals .- 24.0

	Light=Spline1D(xvals,yvals, periodic=true)
	Light2(t)=Light(rem(t,xvals[end]))

	return Light2, df
end

end
