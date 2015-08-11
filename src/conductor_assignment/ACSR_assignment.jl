
""" 2015-08-11. Ported from .m file provided by Jonathon Martin
This function uses a lookup table to find the ACSR (Aluminum Conductor
Steel Reinforced) bare overhead conductor which most closely fits the
current limit of a transmission line.

Input:
`I_lim`: The specified current limit of the transmission line in amps.
`V_base`: The nominal voltage rating of the transmission line in kV.

Outputs:
`D`: Diameter of entire conductor [mm]
`Al_m`: Mass of aluminum portion of conductor [kg/m]
`St_m`: Mass of the steel portion of the conductor [kg/m]
`R`: Resistance of conductor [ohms/m]
`bundle`: Number of conductors per phase

Define ACSR conductor lookup table
Information contained in this table is pulled from
ACSR_datasheet_southwire.pdf which was downloaded from the website of
http://www.southwire.com/ on May 27, 2014

NOTE: If this table is edited, conductors must be sorted by increasing
ampacity ratings in order for this function to work correctly.

The data is organized as follows:
Code Word: naming convention of ACSR conductor
Stranding: aluminum / steel strand count
Diameter: of complete cable [inches] (then converted to [mm])
Aluminum weight: [lbs/1000ft] (then converted to [kg/m])
Steel weight: [lbs/1000ft] (then converted to [kg/m])
Resistance: [ohms/1000ft] (then converted to [ohms/m]) (Assumes an AC current and a conductor temp of 75 degrees C)
Rated Ampacity: [amps] (Assumes a conductor temp of 75 C, ambient temp 25 C, emissivity 0.5, wind 2ft/sec, in sun)
"""
function [D, Al_m, St_m, R, bundle] = ACSR_assignment(I_lim,V_base)

ACSR = [
    #'Turkey'   ,   '6/1',
                          0.198,   24,  12,  0.806,  105;
    #'Swan'     ,   '6/1',
                          0.250,   39,  18,  0.515,  140;
    #'Sparrow'  ,   '6/1',
                          0.316,   62,  29,  0.332,  184;
    #'Robin'    ,   '6/1',
                          0.354,   78,  37,  0.268,  212;
    #'Raven'    ,   '6/1',
                          0.398,   99,  47,  0.217,  242;
    #'Quail'    ,   '6/1',
                          0.447,  124,  59,  0.176,  276;
    #'Pigeon'   ,   '6/1',
                          0.502,  156,  74,  0.144,  315;
    #'Penguin'  ,   '6/1',
                          0.563,  197,  93,  0.119,  357;
    #'Waxwing'  ,  '18/1',
                          0.609,  250,  39, 0.0787,  449;
    #'Partridge',  '26/7',
                          0.642,  251, 115, 0.0779,  475;
    #'Ostrich'  ,  '26/7',
                          0.680,  283, 130, 0.0693,  492;
    #'Merlin'   ,  '18/1',
                          0.684,  315,  49, 0.0625,  519;
    #'Linnet'   ,  '26/7',
                          0.720,  317, 146, 0.0618,  529;
    #'Oriole'   ,  '30/7',
                          0.741,  318, 209, 0.0613,  535;
    #'Chickadee',  '18/1',
                          0.743,  373,  58, 0.0529,  576;
    #'Brant'    ,  '24/7',
                          0.772,  374, 137, 0.0526,  584;
    #'Ibis'     ,  '26/7',
                          0.783,  374, 172, 0.0523,  587;
    #'Lark'     ,  '30/7',
                          0.806,  375, 247, 0.0519,  594;
    #'Pelican'  ,  '18/1',
                          0.814,  447,  70, 0.0442,  646;
    #'Flicker'  ,  '24/7',
                          0.846,  449, 164, 0.0439,  655;
    #'Hen'      ,  '30/7',
                          0.883,  450, 296, 0.0433,  666;
    #'Osprey'   ,  '18/1',
                          0.879,  522,  82, 0.0379,  711;
    #'Parakeet' ,  '24/7',
                          0.914,  524, 192, 0.0376,  721;
    #'Eagle'    ,  '30/7',
                          0.953,  525, 345, 0.0372,  734;
    #'Peacock'  ,  '24/7',
                          0.953,  570, 209, 0.0346,  760;
    #'Wood Duck',  '30/7',
                          0.994,  571, 375, 0.0342,  774;
    #'Rook'     ,  '24/7',
                          0.977,  599, 219, 0.0330,  784;
    #'Scoter'   ,  '30/7',
                          1.019,  600, 395, 0.0325,  798;
    #'Gannet'   ,  '26/7',
                          1.014,  628, 289, 0.0313,  812;
    #'Starling' ,  '26/7',
                          1.051,  674, 310, 0.0292,  849;
    #'Redwing'  , '30/19',
                          1.081,  676, 435, 0.0290,  859;
    #'Coot'     ,  '36/1',
                          1.040,  745,  58, 0.0268,  884;
    #'Drake'    ,  '26/7',
                          1.107,  749, 344, 0.0263,  907;
    #'Mallard'  , '30/19',
                          1.140,  751, 483, 0.0261,  918;
    #'Canary'   ,  '54/7',
                          1.162,  848, 310, 0.0241,  961;
    #'Cardinal' ,  '54/7',
                          1.196,  899, 329, 0.0228,  996;
    #'Curlew'   ,  '54/7',
                          1.245,  973, 356, 0.0211, 1047;
    #'Finch'    , '54/19',
                          1.292, 1053, 375, 0.0197, 1093;
    #'Bunting'  ,  '45/7',
                          1.302, 1123, 219, 0.0182, 1139;
    #'Bittern'  ,  '45/7',
                          1.345, 1198, 234, 0.0171, 1184;
    #'Dipper'   ,  '45/7',
                          1.386, 1273, 248, 0.0162, 1229;
    #'Bobolink' ,  '45/7',
                          1.427, 1348, 263, 0.0153, 1272;
    #'Lapwing'  ,  '45/7',
                          1.504, 1498, 292, 0.0139, 1354;
    #'Chukar'   , '84/19',
                          1.602, 1685, 386, 0.0125, 1453;
    #'Kiwi'     ,  '72/7',
                          1.735, 2051, 249, 0.0106, 1607;
    #'Bluebird' , '84/19',
                          1.762, 2040, 468, 0.0105, 1623;
    ];

# Convert diameter from [inches] to [mm] (1 in = 25.4 mm)
ACSR[:,1] = ACSR[:,1]*25.4

# Convert weights from [lbs/1000ft] to [kg/m] (1 lb = 0.453592 kg and 1000 ft = 304.8 m)
ACSR[:,2:3] = ACSR[:,2:3]*0.453592/304.8

# Convert resistance from [ohms/1000ft] to [ohms/m] (1000ft = 304.8 m)
ACSR[:,4] = ACSR[:,4]/304.8
# Note: resistance values are scaled in order to better match the RTS96
# system data.  These values are tuned so that estimateLineLength.m will
# give a result very close to the lengths provided for the RTS96 system.
resistanceFactor = 2.
ACSR[:,4] = ACSR[:,4]/resistanceFactor

# Store the maximum current rating
maxCurrent = maximum(ACSR[:,5])

# Find conductor which best matches input current limit
# Note that the following search process assumes that the current limits in
# the ACSR table are increasing as the row number increases.

# Determine the most likely conductor bundling configuration based on the
# voltage level and current limit. bundle represents the number of
# conductors per phase.
# This heuristic is based on information in Bergen and Vittal "Power
# Systems Analysis" 2nd ed 2000 p.85, and Weedy and Cory "Electric Power
# Systems" 4th ed 1998 p.129
largeBundle = 1000
if V_base <= 138
    if I_lim <= 1000
        bundle = 1
    elseif I_lim <= maxCurrent*2
        bundle = 2
    else
        # Rating is most likely unrealistic
        bundle = largeBundle
    end
elseif V_base <= 220
    if I_lim <= 1000
        bundle = 1
    elseif I_lim <= maxCurrent*2
        bundle = 2
    elseif I_lim <= maxCurrent*3
        bundle = 3
    else
        # Rating is most likely unrealistic
        bundle = largeBundle
    end
elseif V_base <= 345
    if I_lim <= maxCurrent
        bundle = 1
    elseif I_lim <= maxCurrent*2
        bundle = 2
    elseif I_lim <= maxCurrent*3
        bundle = 3
    else
        # Rating is most likely unrealistic
        bundle = largeBundle
    end
elseif V_base <= 500
    if I_lim <= maxCurrent*2
        bundle = 2
    elseif I_lim <= maxCurrent*4
        bundle = 4
    elseif I_lim <= maxCurrent*6
        bundle = 6
    else
        # Rating is most likely unrealistic
        bundle = largeBundle
    end
else
    if I_lim <= maxCurrent*3
        bundle = 3
    elseif I_lim <= maxCurrent*4
        bundle = 4
    elseif I_lim <= maxCurrent*6
        bundle = 6
    else
        # Rating is most likely unrealistic
        bundle = largeBundle
    end
end

# Initialize search variables
row = 1
numOfRows = length(ACSR[:,1])
found = false

# Perform search for ACSR match
while !found && row <= numOfRows
    if I_lim <= ACSR[row,5]*bundle
        found = true
        row -= 1
        if row < 1
            row = 1
        end
    else
        row += 1
    end
end

# Assign output variables
if !found
    warn("ACSR type not found for line with current rating $I_lim A and voltage rating $V_base kV")
    D = NaN
    Al_m = NaN
    St_m = NaN
    R = NaN
else
    D = ACSR[row,1]
    Al_m = ACSR[row,2]
    St_m = ACSR[row,3]
    R = ACSR[row,4]
end
