# This module is used to achieve a higher variety of different colors for plotting.
# 20 colors are listed in the colors-field of the ColorManager type.
# getColor Method cycles through the array and starts at the begining after reaching the end
# Michael Garstka, 11/24/2016

module colorModule
    export ColorManager, getColor, resetColorCycle

    type ColorManager
        colors::Array{String}
        currentColorIndex::Int64

        function ColorManager()
            colordefs=["#3b44ba",
            "#c8ce24",
            "#a756de",
            "#0ac753",
            "#ff53c4",
            "#01c3ba",
            "#f5218a",
            "#4cd5ff",
            "#f74b3a",
            "#016ed6",
            "#de8500",
            "#8a1e95",
            "#6c4d08",
            "#aca0ff",
            "#a2172f",
            "#abaae1",
            "#ffb38d",
            "#83317e",
            "#9c6f4b",
            "#ff91d1"]
            new(colordefs,1)
        end
    end

    function getColor(CM::colorModule.ColorManager)::String

        totalNumColors = length(CM.colors)

        if CM.currentColorIndex > totalNumColors
            CM.currentColorIndex = 1
        end

        nextColorString = CM.colors[CM.currentColorIndex]
        CM.currentColorIndex = CM.currentColorIndex + 1

        return nextColorString
    end

    function resetColorCycle(CM::colorModule.ColorManager)
        CM.currentColorIndex = 1
        return nothing
    end

end #END MODULE


function formatFloat!(x::Float64,numDecimals::Int64)
    n = 0
    if x < 0.0
        n = 1
    end
    if abs(x) < 10.0
        p = 1
    elseif abs(x) < 100.0 && abs(x) >= 10.0
        p = 2
    elseif abs(x) <  1000.0 && abs(x) >= 100.0
        p = 3
    else
        error("Please enter an absolute value of less than 1000.0!")
    end
    formatX = string(x)
        sum = p+1+numDecimals+n
        len = length(formatX)
        if len < sum
            for i=1:(sum-len)
                formatX = formatX*"0"
            end
        end
    return formatX[1:p+1+numDecimals+n]
end