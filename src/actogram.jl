




module actogram

    using Plots

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
        p=plot(0:48, 0:48, label="", yflip=true) #set up a plot with the right axes
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


        display(p)
        #Now add the circadian predictions





    end



end
