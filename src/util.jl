module util


function unwrap(phaseIn)
    phase=copy(phaseIn)
    for i in 1:length(phase)
        if phase[i]<0.0
            phase[i]+=2π
        end
    end

    for i in 1:length(phase)
        phase[i]+=(i-1)*2π
    end
    return phase
end

end
