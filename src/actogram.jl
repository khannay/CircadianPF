




module actogram

using Plots, Dierckx


    function rectangle(w,x,y)
        Shape(x .+ [0,w,w,0], y .+ [0,0,1,1])
    end


    function getRect(timeon, timeoff, num_days)
        bottom_x=mod(timeon, 24.0)
        bottom_y=num_days-floor(timeon/24.0)-1
        r1=rectangle(timeoff-timeon, bottom_x, bottom_y)
        r2=rectangle(timeoff-timeon, bottom_x+24.0, bottom_y)
        return (r1,r2)
    end


    function makeActogram(tsdf; threshold=10.0)

        num_days=ceil(tsdf.Time[end]/24.0)
        println("Number of days is: $num_days")

        #Add the light schedule to the plot
        p=plot(0:48, 0:48, label="", yflip=true, xtickfontsize=4, ytickfontsize=4, xguidefontsize=6, yguidefontsize=6, titlefontsize=8) #set up a plot with the right axes
        ylims!(p,0,num_days)
        xlims!(p,0,48)
        yticks!(p, 0:num_days)
        xticks!(mod.(0:3:48,24))
        xlabel!("Time of Day")
        ylabel!("Day")
        title!("Actogram")



        lightdata=tsdf.Light_Level
        timedata=tsdf.Time

        lightsOn=false

        if (lightdata[1]>threshold)
            lightsOn=true
            lightStart=timedata[1]
        else
            darkOn=true
            darkStart=timedata[1]
        end

        dayCounter=floor(timedata[1]/24.0)

        for i=1:length(lightdata)
            currentDay=floor(timedata[i]/24.0)

            if (currentDay ≠ dayCounter)
                dayCounter=currentDay
                if (lightsOn)
                    r1,r2=getRect(lightStart, timedata[i], num_days)
                    plot!(p, r1, color=:yellow, label="")
                    plot!(p, r2, color=:yellow, label="")
                    if (i+1<length(timedata))
                        lightStart=timedata[i+1]
                    end
                else
                    r1,r2=getRect(darkStart, timedata[i], num_days)
                    plot!(p, r1, color=:black, label="")
                    plot!(p, r2, color=:black, label="")
                    if (i+1< length(timedata))
                        darkStart=timedata[i+1]
                    end
                end
            end

            if (lightdata[i]<threshold && lightsOn)
                r1,r2=getRect(lightStart, timedata[i-1], num_days)
                plot!(p,r1, color=:yellow, label="")
                plot!(p,r2, color=:yellow, label="")
                lightsOn=false
                darkOn=true
                darkStart=timedata[i]
            end

            if (!(lightsOn) && lightdata[i] >= threshold)
                lightsOn=true
                lightStart=timedata[i]
                darkOn=false
                r1,r2=getRect(darkStart, timedata[i-1], num_days)
                plot!(p,r1, color=:black, label="")
                plot!(p,r2, color=:black, label="")
            end
        end

        dayYvalsDLMO, dlmo_times=getPhaseMarker(tsdf.Time, tsdf.Phase, num_days, marker=1.309)
        dayYvalsCBT, cbt_times=getPhaseMarker(tsdf.Time, tsdf.Phase, num_days, marker=π)

        scatter!(dlmo_times, dayYvalsDLMO, seriestype=:scatter, label="", color=:blue)
        scatter!(dlmo_times .+ 24.0, dayYvalsDLMO, seriestype=:scatter, label="", color=:blue)

        scatter!(cbt_times, dayYvalsCBT, seriestype=:scatter, label="", color=:red, markershape=:cross)
        scatter!(cbt_times .+ 24.0, dayYvalsCBT, seriestype=:scatter, label="", color=:red, markershape=:cross)

        display(p)
        return p
    end

    function getPhaseMarker(Time, Phase, num_days; marker=1.309)
        dlmo_func=Spline1D(Phase, Time)

        real_days=Time[end]/24.0

        if (Phase[1]<marker)
            dlmo_phases=collect(marker:2*π:real_days*2*π)
            dlmo_times=mod.(dlmo_func(dlmo_phases),24.0)
            dayYvalsDLMO=num_days .- collect(0.5:1.0:length(dlmo_times)+0.5)
        else
            dlmo_phases=collect(marker+2*π:2*π:real_days*2*π)
            dlmo_times=mod.(dlmo_func(dlmo_phases),24.0)
            dayYvalsDLMO=num_days .- collect(1.5:1.0:length(dlmo_times)+1.5)
        end

        l1=length(dlmo_times)
        l2=length(dayYvalsDLMO)

        #println("The lengths are $l1 and $l2")

        return dayYvalsDLMO, dlmo_times
    end

end

include("SinglePopModel.jl")
include("HCHSLight.jl")
my_parms=SinglePopModel.sp_parameters()
SinglePopModel.setParameters(my_parms)
L, hchsDataFrame=HCHSLight.readData("./data/hchs-sol-sueno-00579338.csv")
SinglePopModel.setB(L)
last=SinglePopModel.integrateTransients(numdays=200)

last[2]=mod(last[2], 2π)
println("Starting is $last")
tsdf=SinglePopModel.getTS(40*24.0, last)
actogram.makeActogram(tsdf)
